"""
Validate transcriptomic fibre-type annotations against manual IHC labels.
"""

from pathlib import Path

import pandas as pd
import scanpy as sc


# --------------------------------------------------
# Settings
# --------------------------------------------------

annotation_file = Path("manual_annotations.xlsx")
mapping_file = Path("sample_mapping.csv")
anndata_dir = Path("annotated_anndata")
results_dir = Path("validation_results")

prediction_column = "fiber_type_from_myh"
validated_labels = ["type_2a", "type_2b"]

results_dir.mkdir(parents=True, exist_ok=True)


def clean_ids(values):
    """Convert fibre identifiers to comparable strings."""

    return (
        pd.to_numeric(values, errors="coerce")
        .astype("Int64")
        .astype(str)
    )


def load_manual_annotations(sheet_name):
    """Load manual type IIa and IIb labels from one Excel sheet."""

    sheet = pd.read_excel(
        annotation_file,
        sheet_name=sheet_name,
    )

    annotations = pd.DataFrame({
        "cell_id": clean_ids(sheet.iloc[:, 0]),
        "type_2a": sheet["2a"],
        "type_2b": sheet["2b"],
    })

    annotations["manual_label"] = pd.NA
    annotations.loc[
        annotations["type_2a"] == 1,
        "manual_label",
    ] = "type_2a"
    annotations.loc[
        annotations["type_2b"] == 1,
        "manual_label",
    ] = "type_2b"

    return annotations[
        annotations["manual_label"].notna()
    ].copy()


def calculate_metrics(reference, prediction, label):
    """Calculate one-vs-rest precision, recall and F1."""

    tp = ((reference == label) & (prediction == label)).sum()
    fp = ((reference != label) & (prediction == label)).sum()
    fn = ((reference == label) & (prediction != label)).sum()

    precision = tp / (tp + fp) if tp + fp else 0
    recall = tp / (tp + fn) if tp + fn else 0
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 0
    )

    return {
        "tp": int(tp),
        "fp": int(fp),
        "fn": int(fn),
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def main():
    mapping = pd.read_csv(mapping_file)

    required_columns = {"excel_sheet", "anndata_file"}
    missing = required_columns.difference(mapping.columns)

    if missing:
        raise ValueError(
            f"Mapping file is missing columns: {sorted(missing)}"
        )

    sample_results = []
    matched_tables = []
    missing_tables = []

    for row in mapping.itertuples(index=False):
        sample = row.excel_sheet
        h5ad_path = anndata_dir / row.anndata_file

        print(f"Processing {sample}")

        manual = load_manual_annotations(sample)
        adata = sc.read_h5ad(h5ad_path)

        if prediction_column not in adata.obs:
            raise KeyError(
                f"{prediction_column} not found in {h5ad_path.name}"
            )

        adata.obs_names = clean_ids(adata.obs_names)

        predictions = (
            adata.obs[prediction_column]
            .astype(str)
            .str.strip()
            .str.lower()
            .to_dict()
        )

        manual["auto_label"] = manual["cell_id"].map(predictions)

        matched = manual[manual["auto_label"].notna()].copy()
        missing_fibres = manual[manual["auto_label"].isna()].copy()

        matched["sample"] = sample
        missing_fibres["sample"] = sample

        matched_tables.append(matched)
        missing_tables.append(missing_fibres)

        result = {
            "sample": sample,
            "n_manual": len(manual),
            "n_matched": len(matched),
            "n_missing": len(missing_fibres),
            "percent_matched": (
                100 * len(matched) / len(manual)
                if len(manual)
                else 0
            ),
        }

        for label in validated_labels:
            metrics = calculate_metrics(
                matched["manual_label"],
                matched["auto_label"],
                label,
            )

            for metric, value in metrics.items():
                result[f"{label}_{metric}"] = value

        result["macro_f1"] = sum(
            result[f"{label}_f1"]
            for label in validated_labels
        ) / len(validated_labels)

        sample_results.append(result)

    results = pd.DataFrame(sample_results)

    results.to_csv(
        results_dir / "fibre_annotation_metrics_per_sample.csv",
        index=False,
    )

    pd.concat(
        matched_tables,
        ignore_index=True,
    ).to_csv(
        results_dir / "matched_fibres.csv",
        index=False,
    )

    pd.concat(
        missing_tables,
        ignore_index=True,
    ).to_csv(
        results_dir / "missing_fibres.csv",
        index=False,
    )

    summary_rows = []

    for label in validated_labels:
        summary_rows.append({
            "label": label,
            "precision_mean": results[f"{label}_precision"].mean(),
            "precision_sd": results[f"{label}_precision"].std(),
            "recall_mean": results[f"{label}_recall"].mean(),
            "recall_sd": results[f"{label}_recall"].std(),
            "f1_mean": results[f"{label}_f1"].mean(),
            "f1_sd": results[f"{label}_f1"].std(),
            "tp_total": results[f"{label}_tp"].sum(),
            "fp_total": results[f"{label}_fp"].sum(),
            "fn_total": results[f"{label}_fn"].sum(),
        })

    summary = pd.DataFrame(summary_rows)

    summary.to_csv(
        results_dir / "fibre_annotation_metrics_summary.csv",
        index=False,
    )

    pd.DataFrame([{
        "validated_classes": ",".join(validated_labels),
        "macro_f1_mean": results["macro_f1"].mean(),
        "macro_f1_sd": results["macro_f1"].std(),
    }]).to_csv(
        results_dir / "fibre_annotation_macro_f1.csv",
        index=False,
    )

    print(results)
    print(summary)


if __name__ == "__main__":
    main()