import pandas as pd
from scipy.stats import pearsonr, spearmanr

ivanov = pd.read_csv(
    "data/gwas/lifespan_female.pheno",
    sep=r"\s+",
    dtype={"IID": str}
)

ivanov = ivanov[["IID", "S18_1537_F"]].rename(
    columns={"S18_1537_F": "lifespan"}
)

ivanov["lifespan"] = pd.to_numeric(
    ivanov["lifespan"],
    errors="coerce"
)

huang = pd.read_csv(
    "data/gwas/Huangetal_female.25c.mean.pheno",
    sep=r"\s+",
    header=None,
    names=["IID", "lifespan"],
    dtype={"IID": str}
)

huang["lifespan"] = pd.to_numeric(
    huang["lifespan"],
    errors="coerce"
)

ivanov["z"] = (
    ivanov["lifespan"] - ivanov["lifespan"].mean()
) / ivanov["lifespan"].std()

huang["z"] = (
    huang["lifespan"] - huang["lifespan"].mean()
) / huang["lifespan"].std()

combined = pd.concat([
    ivanov[["IID", "z"]],
    huang[["IID", "z"]]
])

combined = (
    combined
    .groupby("IID", as_index=False)["z"]
    .mean()
)

fam = pd.read_csv(
    "data/gwas/tmp/merged_qc.fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "PID", "MID", "SEX", "PHENO"],
    dtype={"IID": str}
)

combined = combined[
    combined["IID"].isin(fam["IID"])
].copy()

combined.insert(0, "FID", "line")
combined = combined[["FID", "IID", "z"]]

combined.to_csv(
    "data/gwas/lifespan_female_combined_z.pheno",
    sep="\t",
    index=False
)

print("Final combined lines:", len(combined))
print()
print(combined.head(10).to_string(index=False))
print()
print(combined["z"].describe())

ivanov_corr = pd.read_csv(
    "data/gwas/lifespan_female.pheno",
    sep=r"\s+",
    dtype={"IID": str}
)

huang_corr = pd.read_csv(
    "data/gwas/Huangetal_female.25c.mean.pheno",
    sep=r"\s+",
    header=None,
    names=["IID", "Huang"],
    dtype={"IID": str}
)

ivanov_corr = ivanov_corr[
    ["IID", "S18_1537_F"]
].rename(
    columns={"S18_1537_F": "Ivanov"}
)

ivanov_corr["Ivanov"] = pd.to_numeric(
    ivanov_corr["Ivanov"],
    errors="coerce"
)

huang_corr["Huang"] = pd.to_numeric(
    huang_corr["Huang"],
    errors="coerce"
)

shared = ivanov_corr.merge(
    huang_corr,
    on="IID",
    how="inner"
)

shared = shared.dropna(
    subset=["Ivanov", "Huang"]
)

pearson_r, pearson_p = pearsonr(
    shared["Ivanov"],
    shared["Huang"]
)

spearman_r, spearman_p = spearmanr(
    shared["Ivanov"],
    shared["Huang"]
)

shared["difference_Huang_minus_Ivanov"] = (
    shared["Huang"] - shared["Ivanov"]
)

print()
print("Shared lines:", len(shared))

print()
print("Correlation:")
print(f"Pearson  r = {pearson_r:.4f}, P = {pearson_p:.4g}")
print(f"Spearman r = {spearman_r:.4f}, P = {spearman_p:.4g}")

print()
print("Difference (Huang - Ivanov):")
print(shared["difference_Huang_minus_Ivanov"].describe())

print()
print("First 20 shared lines:")
print(shared.head(20).to_string(index=False))