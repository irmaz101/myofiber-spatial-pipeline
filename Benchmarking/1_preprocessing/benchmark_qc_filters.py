"""Compare count and density QC decisions, area correlations and spatial patterns."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Settings

input_dir = Path("anndata_mf_density")
output_dir = Path("qc_benchmark")

output_dir.mkdir(parents=True, exist_ok=True)


def analyze_sample(csv_path):
    """Benchmark the two QC filters for one sample."""

    sample = csv_path.stem.replace("_filtering_comparison", "")

    data = pd.read_csv(csv_path, index_col=0)

    required_columns = {
        "total_counts",
        "cell_area",
        "transcript_density",
        "keep_count_filter",
        "keep_density_filter",
        "x_centroid",
        "y_centroid",
    }

    missing = required_columns.difference(data.columns)

    if missing:
        raise ValueError(f"{csv_path.name} is missing columns: {sorted(missing)}")

    keep_count = data["keep_count_filter"].astype(bool)
    keep_density = data["keep_density_filter"].astype(bool)

    removed_count = ~keep_count
    removed_density = ~keep_density
    removed_both = removed_count & removed_density

    n_removed_count = int(removed_count.sum())
    n_removed_density = int(removed_density.sum())
    n_removed_both = int(removed_both.sum())

    summary = {
        "sample": sample,
        "n_objects": len(data),
        "removed_by_count": n_removed_count,
        "removed_by_density": n_removed_density,
        "removed_by_both": n_removed_both,
        "removed_overlap_vs_count": (
            n_removed_both / n_removed_count if n_removed_count else 0
        ),
        "removed_overlap_vs_density": (
            n_removed_both / n_removed_density if n_removed_density else 0
        ),
        "correlation_counts_area": data["total_counts"].corr(data["cell_area"]),
        "correlation_density_area": data["transcript_density"].corr(data["cell_area"]),
        "correlation_counts_density": data["total_counts"].corr(
            data["transcript_density"]
        ),
    }

    data["filter_category"] = "kept_by_both"
    data.loc[removed_count & removed_density, "filter_category"] = "removed_by_both"
    data.loc[removed_count & keep_density, "filter_category"] = "removed_only_by_count"
    data.loc[keep_count & removed_density, "filter_category"] = (
        "removed_only_by_density"
    )

    # Counts versus area
    plt.figure(figsize=(7, 5))
    plt.scatter(data["cell_area"], data["total_counts"], s=5, alpha=0.5)
    plt.xlabel("Object area")
    plt.ylabel("Total transcript count")
    plt.tight_layout()
    plt.savefig(
        output_dir / f"{sample}_counts_vs_area.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    # Density versus area
    plt.figure(figsize=(7, 5))
    plt.scatter(data["cell_area"], data["transcript_density"], s=5, alpha=0.5)
    plt.xlabel("Object area")
    plt.ylabel("Transcript density")
    plt.tight_layout()
    plt.savefig(
        output_dir / f"{sample}_density_vs_area.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    # Spatial distribution of filtering categories
    plt.figure(figsize=(7, 7))

    for category, subset in data.groupby("filter_category"):
        plt.scatter(
            subset["x_centroid"], subset["y_centroid"], s=5, alpha=0.6, label=category
        )

    plt.gca().invert_yaxis()
    plt.xlabel("x centroid")
    plt.ylabel("y centroid")
    plt.legend(markerscale=3, fontsize=8)
    plt.tight_layout()
    plt.savefig(
        output_dir / f"{sample}_filter_categories.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    return summary


def main():
    comparison_files = sorted(input_dir.glob("*_filtering_comparison.csv"))

    if not comparison_files:
        raise FileNotFoundError(
            f"No filtering-comparison files found in " f"{input_dir.resolve()}."
        )

    summaries = []

    for csv_path in comparison_files:
        print(f"Processing {csv_path.name}")

        try:
            summaries.append(analyze_sample(csv_path))

        except Exception as error:
            print(f"Skipped {csv_path.name}: {error}")

    results = pd.DataFrame(summaries)

    if results.empty:
        raise ValueError("No valid filtering-comparison files were processed.")

    results.to_csv(output_dir / "qc_filter_metrics_per_sample.csv", index=False)

    overall = (
        results.select_dtypes("number")
        .agg(["mean", "std", "median"])
        .T.reset_index()
        .rename(columns={"index": "metric"})
    )

    overall.to_csv(output_dir / "qc_filter_metrics_summary.csv", index=False)

    print("\nPer-sample results:")
    print(results)

    print("\nOverall summary:")
    print(overall)


if __name__ == "__main__":
    main()
