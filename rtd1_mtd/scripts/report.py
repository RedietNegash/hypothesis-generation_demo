"""Summarize observed results without cross-dataset mechanism inference."""
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache/matplotlib'))
import json
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def md_table(t):
    # Avoid an optional tabulate dependency.
    rows = ['| '+' | '.join(map(str,t.columns))+' |', '| '+' | '.join(['---']*len(t.columns))+' |']
    rows.extend('| '+' | '.join(map(str,row))+' |' for row in t.itertuples(index=False,name=None))
    return '\n'.join(rows)


def main():
    out = ROOT / 'results'
    bulk = pd.read_csv(out/'bulk/contrast_summary.tsv',sep='\t')
    sensitivity = pd.read_csv(out/'bulk_sensitivity/contrast_summary.tsv',sep='\t')
    comparisons=[]
    for row in bulk.itertuples():
        primary=pd.read_csv(out/f'bulk/{row.contrast}.tsv',sep='\t').set_index('gene_id')
        blocked=pd.read_csv(out/f'bulk_sensitivity/{row.contrast}.tsv',sep='\t').set_index('gene_id')
        shared=primary.index[(primary.padj<.05)&(blocked.reindex(primary.index).padj<.05)]
        comparisons.append({'contrast':row.contrast,'primary_FDR_genes':row.significant,
                            'candidate_block_FDR_genes':int((blocked.padj<.05).sum()),
                            'shared_FDR_genes':len(shared),
                            'shared_same_direction':int((np.sign(primary.loc[shared,'log2FoldChange'])==np.sign(blocked.loc[shared,'log2FoldChange'])).sum())})
    comparison=pd.DataFrame(comparisons)
    comparison.to_csv(out/'bulk_sensitivity/comparison.tsv',sep='\t',index=False)
    qc = pd.read_csv(out/'singlecell/qc_summary.tsv',sep='\t')
    ann = pd.read_csv(ROOT/'data/cluster_annotation.tsv',sep='\t',dtype={'cluster':str})
    cells = pd.read_csv(out/'singlecell/cell_metadata.tsv.gz',sep='\t',index_col=0,dtype={'leiden':str})
    assert set(cells.leiden) == set(ann.cluster)
    labels = ann.set_index('cluster').provisional_cell_type
    cells['provisional_cell_type'] = cells.leiden.map(labels)
    cells.to_csv(out/'singlecell/cell_metadata_annotated.tsv.gz',sep='\t')
    pd.crosstab(cells.provisional_cell_type,cells['sample']).to_csv(out/'singlecell/provisional_celltype_counts.tsv',sep='\t')
    ann.to_csv(out/'singlecell/provisional_annotations.tsv',sep='\t',index=False)
    a = ad.read_h5ad(ROOT/'data/processed/rat_synovium.h5ad',backed='r')
    assert a.obs_names.equals(cells.index)
    coordinates = a.obsm['X_umap']
    fig, ax = plt.subplots(figsize=(10,7))
    for i,label in enumerate(sorted(cells.provisional_cell_type.unique())):
        selected = cells.provisional_cell_type.eq(label).to_numpy()
        ax.scatter(coordinates[selected,0],coordinates[selected,1],s=2,alpha=.65,label=label,color=plt.get_cmap('tab20')(i))
    ax.set(xlabel='UMAP 1',ylabel='UMAP 2',title='Rat synovium: provisional broad cell types\nControl/clodronate liposomes; no peptide treatment')
    ax.legend(loc='upper left',bbox_to_anchor=(1,1),markerscale=4,fontsize=8)
    fig.savefig(out/'singlecell/umap_celltypes.png',dpi=180,bbox_inches='tight');plt.close(fig)
    a.file.close()
    fig,ax=plt.subplots(figsize=(9,4))
    ax.bar(bulk.contrast,bulk.up,color='#287a66',label='Up')
    ax.bar(bulk.contrast,-bulk.down,color='#b75b45',label='Down')
    ax.axhline(0,color='black',linewidth=.6)
    ax.set(ylabel='Genes (down plotted below zero)',title='Human THP-1 RTD-1 response: DESeq2 FDR <0.05')
    ax.tick_params(axis='x',rotation=35); ax.legend()
    fig.tight_layout();fig.savefig(out/'bulk/deg_counts.png',dpi=180);plt.close(fig)
    targets = ['basal_10_vs_0_up','LPS_10_vs_0_down']
    fig,axes=plt.subplots(1,2,figsize=(14,5))
    evidence=[]
    for ax,name in zip(axes,targets):
        t=pd.read_csv(out/f'pathways/{name}.tsv',sep='\t')
        significant=t[t.padj<.05].head(6).iloc[::-1]
        ax.barh(significant.pathway,-np.log10(significant.padj),color='#355c7d')
        ax.set(xlabel='−log10(BH adjusted p)',title=name)
        evidence.append((name,t.iloc[0]))
    fig.tight_layout();fig.savefig(out/'pathways/human_10ug_enrichment.png',dpi=180,bbox_inches='tight');plt.close(fig)
    audit=json.loads((ROOT/'data/processed/rat_2019_audit.json').read_text())
    mapping=pd.read_csv(out/'pathways/mapping_audit.tsv',sep='\t')
    rat=mapping[mapping.analysis.str.startswith('rat_day')]
    text=f'''# RTD-1 / MTD project: Specific Objectives 1 and 2

## Scope and available data

Completed public-data preparation, exploratory gene-expression analysis, and pathway analysis. Network integration, peptide target prediction, AI mechanism ranking, and RTD-1 versus MTD comparison are outside the requested scope. The [dataset inventory](../data/dataset_inventory.tsv) records six dataset/assay entries and their limitations. [Source checksums](../data/source_manifest.json) cover downloaded primary inputs.

The three expression resources remain separate: rat arthritis processed tables, human THP-1 bulk RNA-seq, and rat single-cell synovium. GSE278567 is human, not rat. GSE269120 compares control versus clodronate liposomes in inflamed joints, not healthy versus arthritic joints, and includes no peptide exposure.

## Rat 2019 supplementary results

S2 contains 617 published disease-associated genes; S5 contains 258 day-3 genes and S6 contains 107 day-5 genes. Their union is **{audit['day3_day5_union']} unique genes**, with {audit['day3_day5_overlap']} shared. This does not reconstruct the paper's 340 genes across days 1–5: a day-1 gene-level table is absent from this workbook. Four S5 IDs and one S6 ID are outside S2. These discrepancies are retained in [the audit](../data/processed/rat_2019_audit.json).

The signed fold-change values are preserved verbatim. S5 reports positive Ccl20 and Mmp9 values despite narrative descriptions of suppression; no sign correction or monotonic reversal conclusion was imposed. There are no per-gene adjusted p values or sample-level counts in these three tables. We therefore do not recompute differential expression, manufacture volcano plots, or claim to reproduce the original DESeq2 analysis.

Species-specific Reactome over-representation of S5/S6 used S2 as a **conditional disease-gene background**, excluding out-of-background IDs. Neither day had a pathway with BH FDR <0.05. This restricted test does not reproduce the paper's enrichment against an unavailable assayed transcriptome and does not refute its biological conclusions.

{md_table(rat[['analysis','input_query','annotated_query_in_background','annotated_background','pathways_tested','significant_FDR_005']])}

## Human bulk RNA-seq: GSE278567

24 samples form eight groups (three samples each): 0, 1, 3, or 10 µg/ml RTD-1, with or without LPS, at 20 h. GEO sample descriptions explicitly map count columns to GSM accessions. The original file contains 25,240 genes and fractional Partek E/M estimated counts. We rounded estimated counts to the nearest integer, filtered to ≥10 counts in ≥3 samples, and retained 17,790 genes. This approximation is documented and should be revisited with count-generation details or raw-read reprocessing before confirmatory use. TMM-normalized counts were not used as DESeq2 inputs.

DESeq2 1.52.0 used `~group`, six peptide-versus-vehicle contrasts within LPS strata, and BH FDR 0.05 separately per contrast. Experiment labels alone do not establish matched blocks, so no pairing was assumed. These are unshrunken fold changes, with no extra effect-size cutoff. Independent filtering changes the number of genes eligible for adjusted p values. R package versions are recorded in [sessionInfo](bulk/sessionInfo.txt).

{md_table(bulk)}

These differential-expression counts refer to supplied gene features, retaining separate pseudoautosomal `_PAR_Y` rows. Pathway analysis removes version suffixes and collapses to unique stable Ensembl gene IDs, so foreground/background counts in its audit can be slightly smaller than the DESeq2 feature counts.

PCA separates LPS status on PC1 (73.1% variance) and shows a systematic shift for the third labelled replicate on PC2 (7.4%). This is an unresolved experimental/batch structure: metadata do not confirm a matched block. A sensitivity model `~replicate_label + group` therefore tests adjustment for the candidate block without replacing the primary analysis. The comparison below documents model dependence; pathway results below refer to the primary model only. No samples were removed.

{md_table(comparison)}

[Primary sample PCA](bulk/sample_pca.pdf) and [candidate-block model details](bulk_sensitivity/methods.txt).

The strongest upregulated pathway association at 10 µg/ml without LPS was **{evidence[0][1]['pathway']}** (BH adjusted p = {evidence[0][1]['padj']:.3g}). In the LPS stratum at 10 µg/ml, the strongest downregulated-set association was **{evidence[1][1]['pathway']}** (BH adjusted p = {evidence[1][1]['padj']:.3g}). These are gene-set associations. A neutrophil pathway name in THP-1 cells does not indicate neutrophil abundance, and enrichment alone does not establish pathway activation or inhibition.

![Bulk differential-expression counts](bulk/deg_counts.png)

![Selected human pathway associations](pathways/human_10ug_enrichment.png)

## Rat single-cell RNA-seq: GSE269120

{md_table(qc)}

19,034 cells passed filters: 200–6,000 detected genes, ≤50,000 counts, and <15% mitochondrial counts. Genes detected in ≥3 cells were retained (18,133 genes). Processing used Scanpy 1.11.5, 10,000-count normalization and log1p, 2,000 HVGs selected with sample-aware selection, 40 PCs, an exact 15-neighbor graph over 30 PCs, UMAP, and Leiden resolution 0.6 (seed 0). No condition-removing batch integration was applied.

The 17 clusters were reviewed against canonical marker panels and top expression markers. Provisional broad identities include fibroblasts, endothelial, mural, B, T, monocyte/macrophage, neutrophil, and antigen-presenting myeloid cells. Cluster 7 remains unresolved; cluster 12 has mixed stromal/lymphocyte markers and remains unresolved because residual doublets or ambient RNA could explain the pattern. Clusters 1 and 5 have heat-shock/immediate-early signatures. [Annotation evidence and cautions](singlecell/provisional_annotations.tsv) accompany each cluster; subtype validation remains open.

Cluster-marker pathway enrichment is descriptive: top 50 positive markers with log2FC >0.5 and ≥25% detection, using expressed rat genes as background. Cell-level marker p values do not supply biological replication. There is one sequenced library per condition, so no replicate-based condition differential expression or cell-composition significance tests were performed. This dataset cannot establish RTD-1- or MTD-specific cellular responses. Doublets and ambient RNA were not explicitly corrected.

![Provisional rat cell types](singlecell/umap_celltypes.png)

## Cytokines, proteins, and MTD availability

The 2019 human RA-FLS cytokine assays and 2017 rat arthritis/protease study provide relevant published observations, but no verified machine-readable sample-level cytokine/proteomics matrix was collected from their inspected supplements. The 2021 MTD12813 study and its supplementary PDF were collected; its underlying datasets are available from the corresponding author on reasonable request. No matched MTD transcriptome or peptide-interactor dataset was verified in these resources. No measurements were inferred from figure images, and no author was contacted.

## Statistical and interpretation limits

Reactome release {(ROOT/'data/raw/reactome_version.txt').read_text().strip()} direct Ensembl mappings were used separately for rat and human. No uppercasing-based orthology substitution was used. Tests are one-sided hypergeometric over-representation, pathway size 5–500 after intersecting the assayed background; BH includes all eligible pathways including zero overlaps, separately for each gene list. The [mapping audit](pathways/mapping_audit.tsv) records annotation coverage and excluded IDs. Rat mappings include electronically inferred orthologous pathways, so associations need biological review. Gene-list restriction, incomplete annotation, count rounding, small bulk group sizes, single-cell sampling, and uncorrected technical artifacts limit interpretation. No cross-study or cross-species merging, target ranking, binding inference, or causal mechanism claim was made.

## Primary sources

- [2019 paper](https://doi.org/10.1152/physiolgenomics.00066.2019) and [Figshare supplement](https://doi.org/10.6084/m9.figshare.8855588).
- [Human bulk GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278567).
- [Rat single-cell GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269120) and [single-cell paper](https://doi.org/10.1186/s41232-024-00361-2).
- [2017 arthritis/protein study](https://doi.org/10.1371/journal.pone.0187868).
- [MTD12813 study](https://doi.org/10.1038/s41598-021-02619-y).
- [DESeq2 methods](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), [Scanpy workflow](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html), and [Reactome downloads](https://reactome.org/download-data).
'''
    (out/'analysis_summary.md').write_text(text)
    print('Wrote analysis_summary.md and summary figures')


if __name__ == '__main__':
    main()
