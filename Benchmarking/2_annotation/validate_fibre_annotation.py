"""
Validate transcriptomic fibre-type annotations against manual IHC labels.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support

MANUAL_FILE = Path("manual_annotations.xlsx")
ANNDATA_DIR = Path("annotated_anndata")
OUTPUT_DIR = Path("validation_results")

PREDICTION_COLUMN = "fiber_type_from_myh"
VALIDATED_TYPES = ("type_2a", "type_2b")

VALID_PREDICTIONS = {"type_1", "type_2a", "type_2x", "type_2b", "other"}

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def normalize_id(value):
    """Convert fibre identifiers to comparable strings."""
    value = str(value).strip()
    return value[:-2] if value.endswith(".0") else value


def load_manual_annotations(sheet_name):
    """Load and validate manual IHC annotations."""
    manual = pd.read_excel(MANUAL_FILE, sheet_name=sheet_name)

    manual = manual.rename(columns={manual.columns[0]: "cell_id"})

    required_columns = {"cell_id", "2a", "2b"}
    missing_columns = required_columns.difference(manual.columns)

    if missing_columns:
        raise ValueError(
            f"Missing columns in sheet {sheet_name}: " f"{sorted(missing_columns)}"
        )

    manual = manual.dropna(subset=["cell_id"]).copy()
    manual["cell_id"] = manual["cell_id"].map(normalize_id)

    # Blank IHC fields denote negative labels.
    manual["2a"] = pd.to_numeric(manual["2a"], errors="coerce").fillna(0)

    manual["2b"] = pd.to_numeric(manual["2b"], errors="coerce").fillna(0)

    for column in ["2a", "2b"]:
        invalid_values = set(manual[column].unique()) - {0, 1}

        if invalid_values:
            raise ValueError(
                f"Invalid values in column {column}, "
                f"sheet {sheet_name}: {sorted(invalid_values)}"
            )

    duplicate_ids = manual.loc[
        manual["cell_id"].duplicated(keep=False), "cell_id"
    ].unique()

    if len(duplicate_ids):
        raise ValueError(
            f"Duplicate IDs in sheet {sheet_name}: " f"{duplicate_ids.tolist()}"
        )

    double_positive = manual["2a"].eq(1) & manual["2b"].eq(1)

    if double_positive.any():
        ids = manual.loc[double_positive, "cell_id"].tolist()

        raise ValueError(
            f"Fibres labelled as both IIa and IIb " f"in sheet {sheet_name}: {ids}"
        )

    manual["reference"] = np.select(
        [manual["2a"].eq(1), manual["2b"].eq(1)],
        ["type_2a", "type_2b"],
        default="double_negative",
    )

    return manual[["cell_id", "reference"]]


def load_predictions(h5ad_path):
    """Load transcriptomic fibre-type predictions."""
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Missing AnnData file: {h5ad_path}")

    adata = ad.read_h5ad(h5ad_path)

    if PREDICTION_COLUMN not in adata.obs.columns:
        raise KeyError(
            f"Column '{PREDICTION_COLUMN}' not found " f"in {h5ad_path.name}"
        )

    predictions = pd.DataFrame(
        {
            "cell_id": [normalize_id(value) for value in adata.obs_names],
            "prediction": (
                adata.obs[PREDICTION_COLUMN].astype(str).str.strip().to_numpy()
            ),
        }
    )

    unexpected_labels = set(predictions["prediction"].unique()) - VALID_PREDICTIONS

    if unexpected_labels:
        raise ValueError(
            f"Unexpected fibre-type labels in "
            f"{h5ad_path.name}: "
            f"{sorted(unexpected_labels)}"
        )

    duplicate_ids = predictions.loc[
        predictions["cell_id"].duplicated(keep=False), "cell_id"
    ].unique()

    if len(duplicate_ids):
        raise ValueError(
            f"Duplicate IDs in {h5ad_path.name}: " f"{duplicate_ids.tolist()}"
        )

    return predictions


def calculate_metrics(reference, prediction, fibre_type):
    """Calculate one-vs-rest metrics."""
    y_true = reference.eq(fibre_type)
    y_pred = prediction.eq(fibre_type)

    precision, recall, f1, _ = precision_recall_fscore_support(
        y_true, y_pred, average="binary", zero_division=0
    )

    return {
        "reference_n": int(y_true.sum()),
        "TP": int((y_true & y_pred).sum()),
        "FP": int((~y_true & y_pred).sum()),
        "FN": int((y_true & ~y_pred).sum()),
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def main():
    sample_names = pd.ExcelFile(MANUAL_FILE).sheet_names

    per_sample_results = []
    matched_tables = []

    for sample in sample_names:
        h5ad_path = ANNDATA_DIR / f"{sample}_with_myh_fiber_types.h5ad"

        manual = load_manual_annotations(sample)
        predictions = load_predictions(h5ad_path)

        matched = manual.merge(
            predictions, on="cell_id", how="inner", validate="one_to_one"
        )

        if matched.empty:
            raise ValueError(f"No matching fibre IDs found for " f"sample {sample}.")

        matched["sample"] = sample
        matched_tables.append(matched)

        sample_f1 = []

        # Double-negative fibres are negative observations
        # in both one-vs-rest evaluations.
        for fibre_type in VALIDATED_TYPES:
            metrics = calculate_metrics(
                matched["reference"], matched["prediction"], fibre_type
            )

            per_sample_results.append(
                {"sample": sample, "fibre_type": fibre_type, **metrics}
            )

            sample_f1.append(metrics["f1"])

        macro_f1 = np.mean(sample_f1)

        for result in per_sample_results[-len(VALIDATED_TYPES) :]:
            result["sample_macro_f1"] = macro_f1

    per_sample = pd.DataFrame(per_sample_results)
    matched_all = pd.concat(matched_tables, ignore_index=True)

    class_summary = (
        per_sample.groupby("fibre_type")
        .agg(
            reference_n=("reference_n", "sum"),
            precision_mean=("precision", "mean"),
            precision_sd=("precision", "std"),
            recall_mean=("recall", "mean"),
            recall_sd=("recall", "std"),
            f1_mean=("f1", "mean"),
            f1_sd=("f1", "std"),
        )
        .reset_index()
    )

    macro_per_sample = (
        per_sample[["sample", "sample_macro_f1"]]
        .drop_duplicates()
        .sort_values("sample")
    )

    macro_summary = pd.DataFrame(
        {
            "n_samples": [len(macro_per_sample)],
            "macro_f1_mean": [macro_per_sample["sample_macro_f1"].mean()],
            "macro_f1_sd": [macro_per_sample["sample_macro_f1"].std(ddof=1)],
        }
    )

    confusion = pd.crosstab(
        matched_all["reference"], matched_all["prediction"], margins=True
    )

    reference_counts = (
        matched_all.groupby(["sample", "reference"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    per_sample.to_csv(OUTPUT_DIR / "metrics_per_sample.csv", index=False)

    class_summary.to_csv(OUTPUT_DIR / "metrics_summary.csv", index=False)

    macro_per_sample.to_csv(OUTPUT_DIR / "macro_f1_per_sample.csv", index=False)

    macro_summary.to_csv(OUTPUT_DIR / "macro_f1_summary.csv", index=False)

    confusion.to_csv(OUTPUT_DIR / "confusion_matrix.csv")

    reference_counts.to_csv(OUTPUT_DIR / "reference_counts_per_sample.csv", index=False)

    matched_all.to_csv(OUTPUT_DIR / "matched_fibres.csv", index=False)

    print(class_summary.round(4).to_string(index=False))
    print()
    print(macro_summary.round(4).to_string(index=False))
    print(f"\nResults saved to: " f"{OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
