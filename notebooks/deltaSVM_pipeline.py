import marimo

__generated_with = "0.7.0"
app = marimo.App(width="medium")


@app.cell
def __(mo):
    mo.md(r"""
    # deltaSVM Pipeline
    End-to-end pipeline for scoring SNP allelic effects on TF binding using deltaSVM.
    """)
    return


@app.cell
def __(mo):
    mo.md(r"""## 0. Environment Setup""")
    return


@app.cell
def __():
    import subprocess
    result = subprocess.run(["bash", "-c", """
BASE_DIR="/mnt/hdd_1/rediet/deltaSVM"
mkdir -p $BASE_DIR/scripts

if [ ! -f "$BASE_DIR/scripts/gkmpredict" ]; then
    git clone https://github.com/Dongwon-Lee/lsgkm.git /tmp/lsgkm
    cd /tmp/lsgkm/src && make
    cp /tmp/lsgkm/bin/gkmpredict $BASE_DIR/scripts/
    cp /tmp/lsgkm/bin/gkmtrain $BASE_DIR/scripts/
fi

if [ ! -f "$BASE_DIR/scripts/deltasvm_subset_multi" ]; then
    g++ -O3 -o $BASE_DIR/scripts/deltasvm_subset_multi $BASE_DIR/scripts/deltasvm_subset_multi.cpp
fi

mkdir -p $BASE_DIR/resources/hs38
mkdir -p $BASE_DIR/gkmsvm_models
mkdir -p $BASE_DIR/snp_batches

echo "Setup complete. Make sure the following are in place:"
echo "  - resources/hs38/hs38.fa"
echo "  - resources/models.weights.txt"
echo "  - resources/thresholds.pbs.tsv"
echo "  - resources/thresholds.obs.tsv"
echo "  - gkmsvm_models/*.model.txt"
echo "  - snp_batches/chrN.tsv"
    """], capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)
    return result, subprocess


if __name__ == "__main__":
    app.run()
