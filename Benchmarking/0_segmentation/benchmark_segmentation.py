"""Benchmark predicted masks against a manually annotated reference subset.

Objects are matched one-to-one at the selected IoU threshold. Because the
manual mask may contain only a subset of objects, the script reports detection
recall but not false positives, precision, or F1.

Area ratio is defined as predicted area / manual reference area.
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
    """Load transcript coordinates using the preprocessing coordinate conversion."""
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


def get_objects(mask):
    """Return positive object IDs and their areas."""
    ids, areas = np.unique(mask[mask > 0], return_counts=True)
    return ids.astype(int), dict(zip(ids, areas))


def count_transcripts(mask, transcripts):
    """Count transcripts assigned to each labelled object."""
    valid = (
        transcripts["x_pixel"].between(0, mask.shape[1] - 1)
        & transcripts["y_pixel"].between(0, mask.shape[0] - 1)
    )
    transcripts = transcripts.loc[valid]

    labels = mask[
        transcripts["y_pixel"].to_numpy(),
        transcripts["x_pixel"].to_numpy(),
    ]
    labels = labels[labels > 0]

    return (
        pd.Series(labels)
        .value_counts()
        .rename_axis("object_id")
        .rename("transcript_count")
    )


def calculate_overlaps(manual_mask, predicted_mask):
    """Calculate metrics for all overlapping manual-predicted object pairs."""
    manual_ids, manual_areas = get_objects(manual_mask)
    predicted_ids, predicted_areas = get_objects(predicted_mask)

    manual = manual_mask.ravel().astype(np.int64)
    predicted = predicted_mask.ravel().astype(np.int64)

    overlap = (manual > 0) & (predicted > 0)
    pairs = np.column_stack([manual[overlap], predicted[overlap]])

    columns = [
        "manual_id", "predicted_id", "intersection_pixels",
        "manual_area", "predicted_area", "area_ratio",
        "absolute_area_error", "relative_area_error", "iou", "dice",
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
        manual_id = int(manual_id)
        predicted_id = int(predicted_id)

        area_m = int(manual_areas[manual_id])
        area_p = int(predicted_areas[predicted_id])
        area_error = abs(area_p - area_m)

        records.append({
            "manual_id": manual_id,
            "predicted_id": predicted_id,
            "intersection_pixels": int(intersection),
            "manual_area": area_m,
            "predicted_area": area_p,
            "area_ratio": area_p / area_m,
            "absolute_area_error": area_error,
            "relative_area_error": area_error / area_m,
            "iou": intersection / (area_m + area_p - intersection),
            "dice": 2 * intersection / (area_m + area_p),
        })

    return pd.DataFrame(records), manual_ids, predicted_ids


def match_objects(pair_table, manual_ids, predicted_ids, threshold):
    """Match objects one-to-one, maximizing valid matches and then total IoU."""
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
                manual_pos[row.manual_id],
                predicted_pos[row.predicted_id],
            ] = bonus + row.iou

    matched = {}

    if weights.size:
        rows, columns = linear_sum_assignment(-weights)

        for i, j in zip(rows, columns):
            if weights[i, j] > 0:
                manual_id = int(manual_ids[i])
                predicted_id = int(predicted_ids[j])
                matched[manual_id] = lookup[(manual_id, predicted_id)]

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
                "area_ratio": np.nan,
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
            "area_ratio": float(row.area_ratio),
            "absolute_area_error": int(row.absolute_area_error),
            "relative_area_error": float(row.relative_area_error),
        })

    return pd.DataFrame(records)


def add_transcript_metrics(matches, manual_mask, predicted_mask, transcripts):
    """Calculate transcript-assignment errors for matched objects."""
    detected = matches["status"].eq("detected")

    manual_counts = count_transcripts(manual_mask, transcripts)
    predicted_counts = count_transcripts(predicted_mask, transcripts)

    matches.loc[detected, "manual_transcript_count"] = (
        matches.loc[detected, "manual_id"].map(manual_counts).fillna(0)
    )
    matches.loc[detected, "predicted_transcript_count"] = (
        matches.loc[detected, "predicted_id"].map(predicted_counts).fillna(0)
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
    """Summarize detection and matched-object metrics."""
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
        "iou",
        "dice",
        "area_ratio",
        "absolute_area_error",
        "relative_area_error",
        "absolute_transcript_error",
        "relative_transcript_error",
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

        pairs, manual_ids, predicted_ids = calculate_overlaps(
            manual_mask, predicted_mask
        )

        pairs.to_csv(
            OUTPUT_DIR / f"{model_name}_all_overlapping_pairs.csv",
            index=False,
        )

        matches = match_objects(
            pairs, manual_ids, predicted_ids, IOU_THRESHOLD
        )
        matches = add_transcript_metrics(
            matches, manual_mask, predicted_mask, transcripts
        )

        matches.to_csv(
            OUTPUT_DIR / f"{model_name}_object_metrics.csv",
            index=False,
        )

        summary = summarize_model(
            model_name, matches, len(predicted_ids)
        )
        summaries.append(summary)

        print(
            f"{model_name}: {summary['n_detected']}/"
            f"{summary['n_reference']} detected "
            f"(recall={summary['detection_recall']:.3f})"
        )

    summary_table = pd.DataFrame(summaries)
    summary_table.to_csv(
        OUTPUT_DIR / "segmentation_summary.csv",
        index=False,
    )

    print(f"Results saved in: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
