"""Descriptive rat synovial scRNA-seq analysis; no peptide-response tests."""
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache/matplotlib'))
os.environ.setdefault('NUMBA_CACHE_DIR', str(ROOT / '.cache/numba'))
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PANELS = {
    'macrophage': ['Cd68','Csf1r','Adgre1','Lyz2','C1qa'],
    'fibroblast': ['Col1a1','Col1a2','Dcn','Lum','Pdgfra'],
    'endothelial': ['Pecam1','Cdh5','Kdr','Emcn','Vwf'],
    'T_cell': ['Cd3d','Cd3e','Trac','Lck'],
    'B_cell': ['Cd79a','Cd79b','Ms4a1','Cd19'],
    'neutrophil': ['S100a8','S100a9','Csf3r','Retnlg'],
    'mural': ['Rgs5','Pdgfrb','Acta2','Tagln'],
    'NK_cell': ['Nkg7','Ncr1','Prf1','Gzmb'],
    'mast_cell': ['Cpa3','Kit','Ms4a2','Tpsb2'],
}


def main():
    out = ROOT / 'results/singlecell'
    out.mkdir(parents=True, exist_ok=True)
    arrays = []
    qc_rows = []
    for gsm, name, condition in [('GSM8306635','Cont','control_liposomes'), ('GSM8306636','Clod','clodronate_liposomes')]:
        a = sc.read_10x_h5(ROOT / f'data/raw/{gsm}_{name}_filtered_feature_bc_matrix.h5')
        a.var_names_make_unique()
        a.obs['sample'] = gsm
        a.obs['condition'] = condition
        a.var['mt'] = a.var_names.str.upper().str.startswith('MT-')
        assert a.var['mt'].sum() >= 5, 'Mitochondrial gene annotation missing'
        sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        a.obs.to_csv(out / f'{gsm}_qc_before_filter.tsv', sep='\t')
        keep = ((a.obs.n_genes_by_counts >= 200) & (a.obs.n_genes_by_counts <= 6000)
                & (a.obs.total_counts <= 50000) & (a.obs.pct_counts_mt < 15))
        qc_rows.append({'sample':gsm,'input_cells':a.n_obs,'retained_cells':int(keep.sum()),
                        'mitochondrial_features':int(a.var['mt'].sum())})
        arrays.append(a[keep].copy())
    a = ad.concat(arrays, label='library', keys=['control','clodronate'], index_unique='-', merge='same')
    sc.pp.filter_genes(a, min_cells=3)
    a.layers['counts'] = a.X.copy()
    print(f'QC complete: {a.n_obs} cells, {a.n_vars} genes', flush=True)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    # Keep full expression for marker summaries; reduce to HVGs only for geometry.
    a.raw = a
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor='seurat', batch_key='sample')
    sc.pp.pca(a, n_comps=40, mask_var='highly_variable', svd_solver='arpack', random_state=0)
    print('PCA complete', flush=True)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=0, transformer='sklearn')
    print('Neighbors complete', flush=True)
    sc.tl.umap(a, random_state=0)
    print('UMAP complete', flush=True)
    sc.tl.leiden(a, resolution=0.6, random_state=0, flavor='igraph', n_iterations=2, directed=False)
    sc.tl.rank_genes_groups(a, groupby='leiden', method='wilcoxon', pts=True, use_raw=True)
    markers = sc.get.rank_genes_groups_df(a, group=None)
    # These are cell-level cluster markers, not replicate-based condition DE.
    geneids = a.var['gene_ids'].to_dict()
    markers['ensembl_id'] = markers.names.map(geneids)
    markers.to_csv(out / 'cluster_markers_descriptive.tsv.gz', sep='\t', index=False)
    present = {label:[g for g in genes if g in a.var_names] for label,genes in PANELS.items()}
    cluster_scores = []
    for cluster in a.obs.leiden.cat.categories:
        sel = a.obs.leiden == cluster
        for label, genes in present.items():
            value = float(a[sel, genes].X.mean()) if genes else 0.
            cluster_scores.append({'cluster':cluster,'panel':label,'mean_log1p_expression':value,'genes':';'.join(genes)})
    pd.DataFrame(cluster_scores).to_csv(out / 'cluster_marker_panels.tsv', sep='\t', index=False)
    # Avoid forced cell labels. Marker panels and top markers support manual review.
    pd.crosstab(a.obs.leiden,a.obs['sample']).to_csv(out / 'cluster_sample_counts.tsv', sep='\t')
    pd.DataFrame(qc_rows).to_csv(out / 'qc_summary.tsv', sep='\t', index=False)
    sc.pl.umap(a, color=['sample','leiden'], show=False)
    plt.savefig(out / 'umap.png', dpi=180, bbox_inches='tight'); plt.close('all')
    sc.pl.dotplot(a, present, groupby='leiden', show=False, standard_scale='var')
    plt.savefig(out / 'marker_dotplot.png', dpi=180, bbox_inches='tight'); plt.close('all')
    a.obs.to_csv(out / 'cell_metadata.tsv.gz', sep='\t')
    pd.DataFrame({'symbol':a.var_names,'ensembl_id':a.var.gene_ids.values}).to_csv(out / 'expressed_genes.tsv',sep='\t',index=False)
    a.write_h5ad(ROOT / 'data/processed/rat_synovium.h5ad', compression='gzip')
    (out / 'methods.json').write_text(json.dumps({
        'species':'Rattus norvegicus','peptide_treated':False,
        'samples':2,'sequenced_libraries_per_condition':1,
        'qc':'200–6000 detected genes; <=50000 UMIs; mitochondrial counts <15%; genes in >=3 cells',
        'normalization':'10000 counts per cell, log1p', 'hvg':2000,
        'embedding':'40 PCs; exact sklearn 15 neighbors on 30 PCs; UMAP seed 0; Leiden resolution 0.6 seed 0',
        'annotation':'Unlabeled clusters; marker panels and cluster markers for manual cell-type review.',
        'limitations':'No replicate-based condition DE, no RTD-1/MTD response inference. Residual doublets and ambient RNA not corrected; upper QC bounds only flag potential multiplets.'},indent=2)+'\n')
    print(pd.DataFrame(qc_rows).to_string(index=False))


if __name__ == '__main__':
    main()
