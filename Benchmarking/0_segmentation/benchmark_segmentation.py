"""
Benchmark predicted segmentation masks against manual reference masks.

The script matches each manual object to the overlapping predicted object
with the highest IoU, calculates segmentation and transcript-assignment
metrics, and saves object-level and summary result tables.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from tifffile import imread


# --------------------------------------------------
# Settings
# --------------------------------------------------

manual_mask_file = Path("manual_mask.tif")
predicted_mask_file = Path("predicted_mask.tif")
transcript_file = Path("transcripts.csv")

output_dir = Path("segmentation_benchmark")
iou_threshold = 0.5

output_dir.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# Mask matching
# --------------------------------------------------

def match_objects(manual_mask, predicted_mask):
    """Match each manual object to its highest-IoU predicted object."""

    if manual_mask.shape != predicted_mask.shape:
        raise ValueError(
            f"Mask shapes differ: "
            f"{manual_mask.shape} vs {predicted_mask.shape}"
        )

    manual_flat = manual_mask.ravel().astype(np.int64)
    predicted_flat = predicted_mask.ravel().astype(np.int64)

    manual_ids, manual_sizes = np.unique(
        manual_flat[manual_flat > 0],
        return_counts=True,
    )
    predicted_ids, predicted_sizes = np.unique(
        predicted_flat[predicted_flat > 0],
        return_counts=True,
    )

    manual_size = dict(zip(manual_ids, manual_sizes))
    predicted_size = dict(zip(predicted_ids, predicted_sizes))

    overlap_pixels = (
        (manual_flat > 0)
        & (predicted_flat > 0)
    )

    pairs = np.column_stack([
        manual_flat[overlap_pixels],
        predicted_flat[overlap_pixels],
    ])

    if len(pairs):
        unique_pairs, intersections = np.unique(
            pairs,
            axis=0,
            return_counts=True,
        )
    else:
        unique_pairs = np.empty((0, 2), dtype=int)
        intersections = np.array([], dtype=int)

    best_matches = {}

    for (manual_id, predicted_id), intersection in zip(
        unique_pairs,
        intersections,
    ):
        size_manual = manual_size[manual_id]
        size_predicted = predicted_size[predicted_id]

        union = size_manual + size_predicted - intersection
        iou = intersection / union
        dice = 2 * intersection / (size_manual + size_predicted)

        area_difference = size_predicted - size_manual
        absolute_area_error = abs(area_difference)

        result = {
            "manual_id": int(manual_id),
            "predicted_id": int(predicted_id),
            "intersection_pixels": int(intersection),
            "manual_area": int(size_manual),
            "predicted_area": int(size_predicted),
            "iou": float(iou),
            "dice": float(dice),
            "area_difference": int(area_difference),
            "absolute_area_error": int(absolute_area_error),
            "relative_area_error": float(
                absolute_area_error / size_manual
            ),
            "match_status": "matched",
        }

        if (
            manual_id not in best_matches
            or iou > best_matches[manual_id]["iou"]
        ):
            best_matches[manual_id] = result

    records = []

    for manual_id in manual_ids:
        if manual_id in best_matches:
            records.append(best_matches[manual_id])
        else:
            size_manual = manual_size[manual_id]

            records.append({
                "manual_id": int(manual_id),
                "predicted_id": np.nan,
                "intersection_pixels": 0,
                "manual_area": int(size_manual),
                "predicted_area": 0,
                "iou": 0.0,
                "dice": 0.0,
                "area_difference": -int(size_manual),
                "absolute_area_error": int(size_manual),
                "relative_area_error": 1.0,
                "match_status": "no_overlap",
            })

    matches = pd.DataFrame(records).sort_values("manual_id")
    matches["true_positive"] = matches["iou"] >= iou_threshold

    return matches, manual_ids, predicted_ids, unique_pairs


# --------------------------------------------------
# Transcript assignment
# --------------------------------------------------

def count_transcripts(mask, transcripts):
    """Count transcripts assigned to each non-background mask object."""

    x = pd.to_numeric(
        transcripts["x_location"],
        errors="coerce",
    ).round()

    y = pd.to_numeric(
        transcripts["y_location"],
        errors="coerce",
    ).round()

    valid = (
        x.notna()
        & y.notna()
        & (x >= 0)
        & (x < mask.shape[1])
        & (y >= 0)
        & (y < mask.shape[0])
    )

    x = x.loc[valid].astype(int)
    y = y.loc[valid].astype(int)

    object_ids = mask[
        y.to_numpy(),
        x.to_numpy(),
    ]

    object_ids = object_ids[object_ids > 0]

    return (
        pd.Series(object_ids)
        .value_counts()
        .rename_axis("object_id")
        .rename("transcript_count")
    )


def add_transcript_metrics(
    matches,
    manual_mask,
    predicted_mask,
    transcripts,
):
    """Add transcript counts and assignment errors to matched objects."""

    manual_counts = count_transcripts(
        manual_mask,
        transcripts,
    )

    predicted_counts = count_transcripts(
        predicted_mask,
        transcripts,
    )

    matches["manual_transcript_count"] = (
        matches["manual_id"]
        .map(manual_counts)
        .fillna(0)
        .astype(int)
    )

    matches["predicted_transcript_count"] = (
        matches["predicted_id"]
        .map(predicted_counts)
        .fillna(0)
        .astype(int)
    )

    matches["transcript_count_difference"] = (
        matches["predicted_transcript_count"]
        - matches["manual_transcript_count"]
    )

    matches["absolute_transcript_error"] = (
        matches["transcript_count_difference"].abs()
    )

    matches["relative_transcript_error"] = np.where(
        matches["manual_transcript_count"] > 0,
        matches["absolute_transcript_error"]
        / matches["manual_transcript_count"],
        np.nan,
    )

    return matches


# --------------------------------------------------
# Summary metrics
# --------------------------------------------------

def summarize_results(
    matches,
    manual_ids,
    predicted_ids,
    overlapping_pairs,
):
    """Calculate object-detection and segmentation summary metrics."""

    true_positives = int(matches["true_positive"].sum())
    false_negatives = int(len(manual_ids) - true_positives)

    matched_predicted_ids = set(
        matches.loc[
            matches["true_positive"],
            "predicted_id",
        ]
        .dropna()
        .astype(int)
    )

    false_positives = int(
        len(predicted_ids) - len(matched_predicted_ids)
    )

    precision = (
        true_positives / (true_positives + false_positives)
        if true_positives + false_positives
        else 0
    )

    recall = (
        true_positives / (true_positives + false_negatives)
        if true_positives + false_negatives
        else 0
    )

    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 0
    )

    overlap_table = pd.DataFrame(
        overlapping_pairs,
        columns=["manual_id", "predicted_id"],
    )

    if overlap_table.empty:
        oversegmentation_rate = 0
        undersegmentation_rate = 0
    else:
        oversegmented = (
            overlap_table.groupby("manual_id")["predicted_id"]
            .nunique()
            .gt(1)
            .sum()
        )

        undersegmented = (
            overlap_table.groupby("predicted_id")["manual_id"]
            .nunique()
            .gt(1)
            .sum()
        )

        oversegmentation_rate = (
            oversegmented / len(manual_ids)
            if len(manual_ids)
            else 0
        )

        undersegmentation_rate = (
            undersegmented / len(predicted_ids)
            if len(predicted_ids)
            else 0
        )

    summary = {
        "n_manual_objects": len(manual_ids),
        "n_predicted_objects": len(predicted_ids),
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_negatives": false_negatives,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "iou_threshold": iou_threshold,
        "oversegmentation_rate": oversegmentation_rate,
        "undersegmentation_rate": undersegmentation_rate,
    }

    metric_columns = [
        "iou",
        "dice",
        "absolute_area_error",
        "relative_area_error",
        "absolute_transcript_error",
        "relative_transcript_error",
    ]

    for column in metric_columns:
        summary[f"mean_{column}"] = matches[column].mean()
        summary[f"sd_{column}"] = matches[column].std()
        summary[f"median_{column}"] = matches[column].median()

    return pd.DataFrame([summary])


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    manual_mask = imread(manual_mask_file)
    predicted_mask = imread(predicted_mask_file)
    transcripts = pd.read_csv(transcript_file)

    required_columns = {"x_location", "y_location"}
    missing = required_columns.difference(transcripts.columns)

    if missing:
        raise ValueError(
            f"Transcript file is missing columns: {sorted(missing)}"
        )

    (
        matches,
        manual_ids,
        predicted_ids,
        overlapping_pairs,
    ) = match_objects(
        manual_mask,
        predicted_mask,
    )

    matches = add_transcript_metrics(
        matches,
        manual_mask,
        predicted_mask,
        transcripts,
    )

    summary = summarize_results(
        matches,
        manual_ids,
        predicted_ids,
        overlapping_pairs,
    )

    object_output = (
        output_dir /
        "segmentation_object_metrics.csv"
    )

    summary_output = (
        output_dir /
        "segmentation_summary_metrics.csv"
    )

    matches.to_csv(object_output, index=False)
    summary.to_csv(summary_output, index=False)

    print(f"Saved object-level metrics: {object_output}")
    print(f"Saved summary metrics: {summary_output}")
    print("\nSummary:")
    print(summary.T)


if __name__ == "__main__":
    main()