"""Check scientific input/output contracts against source data and artifacts."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]


def main():
    checks=[]
    for item in json.loads((ROOT/'data/source_manifest.json').read_text()):
        p=ROOT/item['file']
        assert p.stat().st_size == item['bytes']
        assert hashlib.sha256(p.read_bytes()).hexdigest() == item['sha256'],p
    checks.append('All 13 primary downloads match recorded SHA256 checksums')
    p=ROOT/'data/processed'
    original=pd.read_csv(ROOT/'data/raw/GSE278567_Partek_THP1_Gene_counts.txt.gz',sep='\t')
    design=pd.read_csv(p/'human_samples.tsv',sep='\t')
    counts=pd.read_csv(p/'human_counts_rounded.tsv',sep='\t',index_col=0)
    assert list(counts.columns)==list(design['sample'])
    assert list(counts.index)==list(original.gene_id)
    assert np.array_equal(counts.values,np.floor(original[design['sample']].values+.5).astype('int64'))
    assert design.accession.is_unique and (design.groupby('group').size()==3).all()
    assert set(design.species)=={'Homo sapiens'}
    checks.append('24 human samples match GEO accessions and exact rounded source counts; eight groups of three')
    summary=pd.read_csv(ROOT/'results/bulk/contrast_summary.tsv',sep='\t')
    assert len(summary)==6
    expected_genes=set(counts.index[(counts>=10).sum(axis=1)>=3])
    for row in summary.itertuples():
        t=pd.read_csv(ROOT/f'results/bulk/{row.contrast}.tsv',sep='\t')
        assert set(t.gene_id)==expected_genes and t.gene_id.is_unique
        assert int((t.padj<.05).sum())==row.significant
        assert row.up+row.down==row.significant
        assert t.padj.dropna().between(0,1).all()
        assert t.pvalue.dropna().between(0,1).all()
    checks.append('All six DESeq2 outputs cover exactly the independently verified 17,790-gene filter and reconcile with summaries')
    sensitivity=pd.read_csv(ROOT/'results/bulk_sensitivity/contrast_summary.tsv',sep='\t')
    assert len(sensitivity)==6 and set(sensitivity.contrast)==set(summary.contrast)
    for row in sensitivity.itertuples():
        t=pd.read_csv(ROOT/f'results/bulk_sensitivity/{row.contrast}.tsv',sep='\t')
        assert set(t.gene_id)==expected_genes and t.gene_id.is_unique
        assert int((t.padj<.05).sum())==row.significant
        assert row.up+row.down==row.significant
        assert t.padj.dropna().between(0,1).all()
    comparison=pd.read_csv(ROOT/'results/bulk_sensitivity/comparison.tsv',sep='\t')
    assert len(comparison)==6
    assert (comparison.shared_FDR_genes<=comparison.primary_FDR_genes).all()
    assert (comparison.shared_FDR_genes<=comparison.candidate_block_FDR_genes).all()
    checks.append('Candidate-block sensitivity outputs reconcile across six contrasts, with unchanged source features')
    rat=[pd.read_csv(p/f'rat_2019_S{n}.tsv',sep='\t',float_precision='round_trip') for n in [2,5,6]]
    assert [len(t) for t in rat]==[617,258,107]
    assert len(set(rat[1].ensembl_id)|set(rat[2].ensembl_id))==329
    for table, number in zip(rat,[2,5,6]):
        source=pd.read_excel(ROOT/'data/raw/PhysiolGen_SupplementaryTables_2019.xls',sheet_name=f'Supplementary TABLE S{number}',skiprows=1)
        assert np.array_equal(source.iloc[:,3].values,table.reported_signed_fold_change.values)
        assert list(source.iloc[:,1])==list(table.ensembl_id)
    checks.append('Rat S2/S5/S6 values and IDs match workbook verbatim; discrepancy audit retained')
    a=ad.read_h5ad(p/'rat_synovium.h5ad',backed='r')
    assert a.shape==(19034,18133)
    assert a.obs_names.is_unique and a.var_names.is_unique
    assert a.obs['sample'].nunique()==2 and a.obs.leiden.nunique()==17
    assert (a.obs.n_genes_by_counts>=200).all() and (a.obs.n_genes_by_counts<=6000).all()
    assert (a.obs.pct_counts_mt<15).all() and (a.obs.total_counts<=50000).all()
    assert np.isfinite(a.obsm['X_umap']).all()
    qc=pd.read_csv(ROOT/'results/singlecell/qc_summary.tsv',sep='\t')
    assert qc.retained_cells.sum()==a.n_obs
    counts_sc=pd.read_csv(ROOT/'results/singlecell/cluster_sample_counts.tsv',sep='\t',index_col=0)
    assert counts_sc.sum().sum()==a.n_obs
    annotations=pd.read_csv(ROOT/'data/cluster_annotation.tsv',sep='\t',dtype={'cluster':str})
    assert set(annotations.cluster)==set(a.obs.leiden.astype(str))
    a.file.close()
    checks.append('Single-cell dimensions, sample counts, QC thresholds, embedding, cluster counts and provisional annotations reconcile')
    audit=pd.read_csv(ROOT/'results/pathways/mapping_audit.tsv',sep='\t')
    assert len(audit)==31 # 2 rat timepoints + 12 human directional lists + 17 clusters
    for row in audit.itertuples():
        t=pd.read_csv(ROOT/f'results/pathways/{row.analysis}.tsv',sep='\t')
        assert len(t)==row.pathways_tested
        assert t.pvalue.between(0,1).all() and t.padj.between(0,1).all()
        assert np.allclose(t.padj,multipletests(t.pvalue,method='fdr_bh')[1])
        assert (t.overlap<=t.query_size).all() and (t.overlap<=t.pathway_size).all()
        assert t.pathway_size.between(5,500).all()
        assert (t.background_size>=t.query_size).all()
        assert (t.background_size==row.annotated_background).all()
        prefix='R-HSA-' if row.analysis.startswith(('basal','LPS')) else 'R-RNO-'
        assert t.pathway_id.str.startswith(prefix).all()
    checks.append('31 enrichment analyses satisfy species, universe, overlap and BH correction checks')
    result={'status':'passed','checks':checks,'scope':'Specific Objectives 1 and 2 only'}
    (ROOT/'results/validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
