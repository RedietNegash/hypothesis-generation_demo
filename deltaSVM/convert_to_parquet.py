import pandas as pd

df = pd.read_csv("out/summary.pred.tsv", sep="\t")


df = df.rename(columns={
    "snp": "variant",
    "tf": "TfName",
    "preferred_allele": "Effect",
    "deltaSVM_score": "Score"
})

df["rsId"] = "."
df["TfId"] = "."


df = df[df["Effect"].isin(["Gain", "Loss"])]
df = df[["variant", "rsId", "TfId", "TfName", "Effect", "Score"]]

print(df)
print(f"\nTotal rows after filtering: {len(df)}")


df.to_parquet("out/deltasvm_results.parquet", index=False)
print("Saved to out/deltasvm_results.parquet")