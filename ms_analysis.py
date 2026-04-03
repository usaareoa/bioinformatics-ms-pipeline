import pandas as pd
import os
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs("outputs", exist_ok=True)
dataset_folder = "datasets"
mapping = pd.read_csv("genome_mapping.txt", sep=",")
mapping = mapping.rename(columns={"Gene stable ID": "gene", "HGNC symbol": "gene_name"})

targets = ["S100A8", "S100A9", "ALOX15B"]

datasets = [
    "GSE126802", "GSE38010", "GSE180759a", "GSE180759b", "GSE32915", "GSE100297",
    "GSE211358", "GSE126427", "GSE244313"
]

tissue_map = {
    "GSE126802": "CNS", "GSE38010": "CNS", "GSE180759a": "CNS",
    "GSE180759b": "CNS", "GSE32915": "CNS", "GSE100297": "CNS",
    "GSE211358": "blood", "GSE126427": "blood", "GSE244313": "blood"
}

def is_log_transformed(expr_df):
    sample_vals = expr_df.select_dtypes(include=[np.number]).values.flatten()
    sample_vals = sample_vals[~np.isnan(sample_vals)]
    sample_vals = np.random.choice(sample_vals, min(1000, len(sample_vals)), replace=False)
    return np.median(sample_vals) < 50

results = []

for gse in datasets:
    filepath_expr = os.path.join(dataset_folder, f"{gse}_expression.csv")
    filepath_meta = os.path.join(dataset_folder, f"{gse}_metadata.csv")

    if not os.path.exists(filepath_expr) or not os.path.exists(filepath_meta):
        print(f"Skipping {gse} — file not found")
        continue

    expr = pd.read_csv(filepath_expr)
    meta = pd.read_csv(filepath_meta)

    already_log = is_log_transformed(expr)
    print(f"{gse} — log-transformed: {already_log}")

    expr = expr.merge(mapping, on="gene", how="left")
    expr = expr.set_index("gene_name")

    ms_samples  = meta[meta["disease"] == "multiple sclerosis"]["sample_name"].tolist()
    ctl_samples = meta[meta["disease"] == "healthy"]["sample_name"].tolist()

    ms_expr  = expr[[s for s in ms_samples  if s in expr.columns]]
    ctl_expr = expr[[s for s in ctl_samples if s in expr.columns]]

    for gene in targets:
        if gene in expr.index:
            ms_vals  = ms_expr.loc[gene].dropna()
            ctl_vals = ctl_expr.loc[gene].dropna()

            ms_mean  = ms_vals.mean()
            ctl_mean = ctl_vals.mean()

            if already_log:
                log2fc = ms_mean - ctl_mean
            else:
                log2fc = np.log2(ms_mean / ctl_mean) if ctl_mean > 0 else 0

            if len(ms_vals) > 1 and len(ctl_vals) > 1:
                t_stat, p_val = stats.ttest_ind(ms_vals, ctl_vals)
            else:
                t_stat, p_val = np.nan, np.nan

            results.append({
                "dataset":  gse,
                "tissue":   tissue_map[gse],
                "gene":     gene,
                "ms_mean":  round(ms_mean, 4),
                "ctl_mean": round(ctl_mean, 4),
                "log2FC":   round(log2fc, 4),
                "ms_n":     len(ms_vals),
                "ctl_n":    len(ctl_vals),
                "p_value":  round(p_val, 4) if not np.isnan(p_val) else "NA"
            })
        else:
            print(f"{gene} not found in {gse}")

results_df = pd.DataFrame(results)
results_df["log2FC"] = pd.to_numeric(results_df["log2FC"], errors="coerce")

print("\n--- Per-dataset results ---")
print(results_df.to_string())

summary_rows = []
for gene in targets:
    for tissue in ["CNS", "blood", "all"]:
        if tissue == "all":
            gene_df = results_df[results_df["gene"] == gene].copy()
        else:
            gene_df = results_df[(results_df["gene"] == gene) & (results_df["tissue"] == tissue)].copy()

        if len(gene_df) == 0:
            continue

        gene_df["p_value"] = pd.to_numeric(gene_df["p_value"], errors="coerce")

        summary_rows.append({
            "gene":           gene,
            "tissue":         tissue,
            "datasets_total": len(gene_df),
            "upregulated":    (gene_df["log2FC"] > 0).sum(),
            "downregulated":  (gene_df["log2FC"] < 0).sum(),
            "mean_log2FC":    round(gene_df["log2FC"].mean(), 4),
            "median_log2FC":  round(gene_df["log2FC"].median(), 4),
            "p<0.05_count":   (gene_df["p_value"] < 0.05).sum()
        })

summary_df = pd.DataFrame(summary_rows)
print("\n--- Summary by gene and tissue ---")
print(summary_df.to_string())

results_df.to_csv("outputs/results.csv", index=False)
summary_df.to_csv("outputs/summary.csv", index=False)
print("\nSaved to results.csv and summary.csv")



sns.set_theme(style="whitegrid", font_scale=1.1)
palette = {"CNS": "#2E86AB", "blood": "#E84855"}

#  log2FC per dataset per gene colored for each tissue
fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=False)
fig.suptitle("log2 Fold Change (MS vs Control) Across Datasets", fontsize=14, fontweight="bold", y=1.02)

for ax, gene in zip(axes, targets):
    gene_df = results_df[results_df["gene"] == gene].copy()
    colors  = [palette[t] for t in gene_df["tissue"]]

    bars = ax.bar(gene_df["dataset"], gene_df["log2FC"], color=colors, edgecolor="white", linewidth=0.5)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_title(gene, fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("log2FC" if gene == targets[0] else "")
    ax.tick_params(axis="x", rotation=45)


    for bar, tissue in zip(bars, gene_df["tissue"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.005,
            tissue[0].upper(),
            ha="center", va="bottom", fontsize=8, color="gray"
        )


cns_patch   = mpatches.Patch(color=palette["CNS"],   label="CNS")
blood_patch = mpatches.Patch(color=palette["blood"], label="Blood")
fig.legend(handles=[cns_patch, blood_patch], loc="upper right", bbox_to_anchor=(1.02, 1.0))

plt.tight_layout()
plt.savefig("outputs/plot1_log2fc_per_dataset.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved plot1_log2fc_per_dataset.png")


summary_plot = summary_df[summary_df["tissue"] != "all"].copy()

fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
fig.suptitle("Datasets Showing Upregulation vs Downregulation by Tissue", fontsize=13, fontweight="bold")

x      = np.arange(2)
width  = 0.35

for ax, gene in zip(axes, targets):
    gene_data = summary_plot[summary_plot["gene"] == gene]
    cns_row   = gene_data[gene_data["tissue"] == "CNS"].iloc[0]
    blood_row = gene_data[gene_data["tissue"] == "blood"].iloc[0]

    up_vals   = [cns_row["upregulated"],   blood_row["upregulated"]]
    down_vals = [cns_row["downregulated"], blood_row["downregulated"]]

    ax.bar(x - width/2, up_vals,   width, label="Upregulated",   color="#2ECC71", edgecolor="white")
    ax.bar(x + width/2, down_vals, width, label="Downregulated", color="#E74C3C", edgecolor="white")

    ax.set_title(gene, fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(["CNS", "Blood"])
    ax.set_ylabel("Number of datasets" if gene == targets[0] else "")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

axes[0].legend(fontsize=9)
plt.tight_layout()
plt.savefig("outputs/plot2_updown_summary.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved plot2_updown_summary.png")