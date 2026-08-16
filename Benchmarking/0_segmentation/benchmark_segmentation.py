"""Benchmark predicted masks against a manually annotated reference subset.

Objects are matched one-to-one at the selected IoU threshold. Because the
manual mask may contain only a subset of objects, the script reports detection
recall but not false positives, precision, or F1.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from tifffile import imread


# Settings
MANUAL_MASK = Path("manual_mask.tif")
TRANSCRIPT_FILE = Path("transcripts.csv")
PREDICTED_MASKS = {
    "default": Path("default_mask.tif"),
    "custom": Path("custom_mask.tif"),
}
IOU_THRESHOLD = 0.5
OUTPUT_DIR = Path("segmentation_benchmark")


def load_transcripts(path):
    """Load coordinates using the same floor conversion as preprocessing."""
    transcripts = pd.read_csv(path)
    required = {"x_location", "y_location"}
    missing = required.difference(transcripts.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")

    x = pd.to_numeric(transcripts["x_location"], errors="coerce")
    y = pd.to_numeric(transcripts["y_location"], errors="coerce")
    valid = np.isfinite(x) & np.isfinite(y)

    transcripts = transcripts.loc[valid].copy()
    transcripts["x_pixel"] = np.floor(x.loc[valid]).astype(int)
    transcripts["y_pixel"] = np.floor(y.loc[valid]).astype(int)
    return transcripts


def count_transcripts(mask, transcripts):
    """Count transcripts assigned to each labelled object."""
    valid = (
        transcripts["x_pixel"].between(0, mask.shape[1] - 1)
        & transcripts["y_pixel"].between(0, mask.shape[0] - 1)
    )
    inside = transcripts.loc[valid]
    object_ids = mask[
        inside["y_pixel"].to_numpy(),
        inside["x_pixel"].to_numpy(),
    ]
    return (
        pd.Series(object_ids[object_ids > 0])
        .value_counts()
        .rename_axis("object_id")
        .rename("transcript_count")
    )


def calculate_pairwise_overlaps(manual_mask, predicted_mask):
    """Calculate metrics for every overlapping manual-predicted pair."""
    manual = manual_mask.ravel().astype(np.int64)
    predicted = predicted_mask.ravel().astype(np.int64)

    manual_ids, manual_areas = np.unique(
        manual[manual > 0], return_counts=True
    )
    predicted_ids, predicted_areas = np.unique(
        predicted[predicted > 0], return_counts=True
    )
    manual_area = dict(zip(manual_ids, manual_areas))
    predicted_area = dict(zip(predicted_ids, predicted_areas))

    overlap = (manual > 0) & (predicted > 0)
    pairs = np.column_stack([manual[overlap], predicted[overlap]])
    columns = [
        "manual_id", "predicted_id", "intersection_pixels",
        "manual_area", "predicted_area", "iou", "dice",
        "absolute_area_error", "relative_area_error",
    ]
    if not len(pairs):
        return pd.DataFrame(columns=columns), manual_ids, predicted_ids

    unique_pairs, intersections = np.unique(
        pairs, axis=0, return_counts=True
    )
    records = []
    for (manual_id, predicted_id), intersection in zip(
        unique_pairs, intersections
    ):
        area_m = manual_area[manual_id]
        area_p = predicted_area[predicted_id]
        area_error = abs(area_p - area_m)
        records.append({
            "manual_id": int(manual_id),
            "predicted_id": int(predicted_id),
            "intersection_pixels": int(intersection),
            "manual_area": int(area_m),
            "predicted_area": int(area_p),
            "iou": intersection / (area_m + area_p - intersection),
            "dice": 2 * intersection / (area_m + area_p),
            "absolute_area_error": int(area_error),
            "relative_area_error": area_error / area_m,
        })
    return pd.DataFrame(records), manual_ids, predicted_ids


def match_objects_one_to_one(pair_table, manual_ids, predicted_ids, threshold):
    """Maximize valid match count, then total IoU, with unique assignments."""
    manual_pos = {object_id: i for i, object_id in enumerate(manual_ids)}
    predicted_pos = {
        object_id: i for i, object_id in enumerate(predicted_ids)
    }
    weights = np.zeros((len(manual_ids), len(predicted_ids)))
    bonus = min(len(manual_ids), len(predicted_ids)) + 1
    lookup = {}

    for row in pair_table.itertuples(index=False):
        key = (int(row.manual_id), int(row.predicted_id))
        lookup[key] = row
        if row.iou >= threshold:
            weights[
                manual_pos[row.manual_id], predicted_pos[row.predicted_id]
            ] = bonus + row.iou

    matched = {}
    if weights.size:
        rows, columns = linear_sum_assignment(-weights)
        for row_index, column_index in zip(rows, columns):
            if weights[row_index, column_index] > 0:
                key = (
                    int(manual_ids[row_index]),
                    int(predicted_ids[column_index]),
                )
                matched[key[0]] = lookup[key]

    records = []
    for manual_id in map(int, manual_ids):
        if manual_id not in matched:
            records.append({
                "manual_id": manual_id,
                "predicted_id": np.nan,
                "status": "missed",
                "iou": np.nan,
                "dice": np.nan,
                "manual_area": np.nan,
                "predicted_area": np.nan,
                "absolute_area_error": np.nan,
                "relative_area_error": np.nan,
            })
            continue

        row = matched[manual_id]
        records.append({
            "manual_id": manual_id,
            "predicted_id": int(row.predicted_id),
            "status": "detected",
            "iou": float(row.iou),
            "dice": float(row.dice),
            "manual_area": int(row.manual_area),
            "predicted_area": int(row.predicted_area),
            "absolute_area_error": int(row.absolute_area_error),
            "relative_area_error": float(row.relative_area_error),
        })
    return pd.DataFrame(records)


def add_transcript_metrics(matches, manual_mask, predicted_mask, transcripts):
    """Add transcript-assignment errors to detected pairs."""
    detected = matches["status"].eq("detected")
    manual_counts = count_transcripts(manual_mask, transcripts)
    predicted_counts = count_transcripts(predicted_mask, transcripts)

    for output, identifier, counts in [
        ("manual_transcript_count", "manual_id", manual_counts),
        ("predicted_transcript_count", "predicted_id", predicted_counts),
    ]:
        matches.loc[detected, output] = (
            matches.loc[detected, identifier].map(counts).fillna(0).astype(int)
        )

    difference = (
        matches.loc[detected, "predicted_transcript_count"]
        - matches.loc[detected, "manual_transcript_count"]
    )
    matches.loc[detected, "transcript_count_difference"] = difference
    matches.loc[detected, "absolute_transcript_error"] = difference.abs()

    valid = detected & matches["manual_transcript_count"].gt(0)
    matches.loc[valid, "relative_transcript_error"] = (
        matches.loc[valid, "absolute_transcript_error"]
        / matches.loc[valid, "manual_transcript_count"]
    )
    return matches


def summarize_model(model_name, matches, n_predicted):
    """Create one summary row from valid one-to-one matches."""
    detected = matches[matches["status"].eq("detected")]
    n_reference = len(matches)
    summary = {
        "model": model_name,
        "n_reference": n_reference,
        "n_predicted_total": n_predicted,
        "n_detected": len(detected),
        "n_missed": n_reference - len(detected),
        "detection_recall": len(detected) / n_reference
        if n_reference else np.nan,
    }
    metrics = [
        "iou", "dice", "absolute_area_error", "relative_area_error",
        "absolute_transcript_error", "relative_transcript_error",
    ]
    for metric in metrics:
        summary.update({
            f"mean_{metric}": detected[metric].mean(),
            f"sd_{metric}": detected[metric].std(ddof=1),
            f"median_{metric}": detected[metric].median(),
        })
    return summary


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manual_mask = imread(MANUAL_MASK)
    transcripts = load_transcripts(TRANSCRIPT_FILE)
    summaries = []

    for model_name, mask_path in PREDICTED_MASKS.items():
        predicted_mask = imread(mask_path)
        if manual_mask.shape != predicted_mask.shape:
            raise ValueError(
                f"Mask shapes differ for {model_name}: "
                f"{manual_mask.shape} vs {predicted_mask.shape}"
            )

        pairs, manual_ids, predicted_ids = calculate_pairwise_overlaps(
            manual_mask, predicted_mask
        )
        pairs.to_csv(
            OUTPUT_DIR / f"{model_name}_all_overlapping_pairs.csv", index=False
        )

        matches = match_objects_one_to_one(
            pairs, manual_ids, predicted_ids, IOU_THRESHOLD
        )
        matches = add_transcript_metrics(
            matches, manual_mask, predicted_mask, transcripts
        )
        matches.to_csv(
            OUTPUT_DIR / f"{model_name}_object_metrics.csv", index=False
        )

        summary = summarize_model(model_name, matches, len(predicted_ids))
        summaries.append(summary)
        print(
            f"{model_name}: {summary['n_detected']}/{summary['n_reference']} "
            f"detected (recall={summary['detection_recall']:.3f})"
        )

    summary_table = pd.DataFrame(summaries)
    summary_table.to_csv(
        OUTPUT_DIR / "segmentation_summary.csv", index=False
    )
    print(f"Results saved in: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
