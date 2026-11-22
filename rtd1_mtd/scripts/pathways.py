"""Species-specific Reactome over-representation, restricted to Objective 2."""
from pathlib import Path
from collections import defaultdict
import json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/pathways'
OUT.mkdir(parents=True, exist_ok=True)


def read_mappings():
    maps = defaultdict(lambda: defaultdict(set))
    names = {}
    records = []
    with (ROOT / 'data/raw/Ensembl2Reactome.txt').open() as stream:
        for line in stream:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 6:
                raise ValueError('Incomplete Reactome mapping line')
            gene, pathway, url, name, evidence, species = fields
            if species in ['Homo sapiens', 'Rattus norvegicus'] and gene.startswith(('ENSG','ENSRNOG')):
                maps[species][pathway].add(gene.split('.')[0])
                names[pathway] = name
                records.append((gene,pathway,name,evidence,species))
    pd.DataFrame(records,columns=['ensembl_id','pathway_id','pathway','evidence','species']).drop_duplicates().to_csv(ROOT / 'data/processed/reactome_gene_pathways.tsv.gz',sep='\t',index=False)
    return maps, names


def ora(query, background, pathways, names, label, min_size=5):
    # Explicitly restrict both query and pathways to the same assayed universe.
    query_raw = {str(x).split('.')[0] for x in query if pd.notna(x)}
    background_raw = {str(x).split('.')[0] for x in background if pd.notna(x)}
    annotated = set().union(*pathways.values())
    universe = background_raw & annotated
    foreground = query_raw & universe
    audit = {'analysis':label,'input_query':len(query_raw),'input_background':len(background_raw),
             'annotated_background':len(universe),'annotated_query_in_background':len(foreground),
             'query_outside_background':len(query_raw-background_raw),
             'query_missing_annotation':len((query_raw&background_raw)-annotated)}
    results = []
    for pid, genes in pathways.items():
        members = genes & universe
        if not min_size <= len(members) <= 500:
            continue
        overlap = foreground & members
        p = float(hypergeom.sf(len(overlap)-1, len(universe), len(members), len(foreground))) if foreground else 1.
        results.append({'pathway_id':pid,'pathway':names[pid], 'overlap':len(overlap),
                        'query_size':len(foreground), 'pathway_size':len(members),
                        'background_size':len(universe),'pvalue':p,
                        'genes':';'.join(sorted(overlap))})
    result = pd.DataFrame(results)
    assert len(result), f'No testable pathways for {label}'
    result['padj'] = multipletests(result.pvalue, method='fdr_bh')[1]
    # All size-eligible pathways, including zero overlaps, enter BH correction.
    result.sort_values(['padj','pvalue']).to_csv(OUT / f'{label}.tsv',sep='\t',index=False)
    audit['pathways_tested'] = len(result)
    audit['significant_FDR_005'] = int((result.padj < .05).sum())
    return audit


def main():
    maps, names = read_mappings()
    audits = []
    p = ROOT / 'data/processed'
    rat_bg = pd.read_csv(p / 'rat_2019_S2.tsv',sep='\t').ensembl_id
    for day, table in [(3,5),(5,6)]:
        rat = pd.read_csv(p / f'rat_2019_S{table}.tsv',sep='\t')
        # Conditional test within the published disease signature only. No
        # whole-transcriptome background was released in the supplementary XLS.
        audits.append(ora(rat.ensembl_id,rat_bg,maps['Rattus norvegicus'],names,f'rat_day{day}_conditional'))
    files = sorted((ROOT / 'results/bulk').glob('*_vs_0.tsv'))
    if len(files) != 6:
        raise ValueError('Run all six DESeq2 contrasts before pathway analysis')
    for file in files:
        t = pd.read_csv(file,sep='\t')
        # Independent filtering defines which genes could receive adjusted p.
        background = t.loc[t.padj.notna(),'gene_id']
        for direction, sign in [('up',1),('down',-1)]:
            genes = t.loc[(t.padj<.05)&(t.log2FoldChange*sign>0),'gene_id']
            audits.append(ora(genes,background,maps['Homo sapiens'],names,file.stem+'_'+direction))
    sc = ROOT / 'results/singlecell'
    markers = pd.read_csv(sc / 'cluster_markers_descriptive.tsv.gz',sep='\t',dtype={'group':str})
    background = pd.read_csv(sc / 'expressed_genes.tsv',sep='\t').ensembl_id
    for cluster, frame in markers.groupby('group'):
        # Top 50 positive cluster markers with >=25% within-cluster detection.
        # No condition effect is tested, and no cell-level p value is used as
        # evidence for a biological treatment effect.
        selected = frame[(frame.logfoldchanges>.5)&(frame.pct_nz_group>=.25)].sort_values('scores',ascending=False).head(50)
        audits.append(ora(selected.ensembl_id,background,maps['Rattus norvegicus'],names,f'rat_cluster_{cluster}_descriptive'))
    pd.DataFrame(audits).to_csv(OUT / 'mapping_audit.tsv',sep='\t',index=False)
    (OUT / 'methods.json').write_text(json.dumps({
        'database':'Reactome direct Ensembl mappings; species-specific; no inferred human-symbol matching',
        'release':(ROOT/'data/raw/reactome_version.txt').read_text().strip(),
        'identifier_handling':'Remove Ensembl version suffix and collapse to unique stable gene IDs; primary and _PAR_Y features map to the same gene for enrichment only. DESeq2 retains the distinct source features.',
        'test':'one-sided hypergeometric over-representation, BH across every size-eligible pathway per gene list',
        'pathway_size':'5–500 genes after universe intersection; all overlaps including zero tested',
        'rat_2019':'conditional enrichment among S2 disease genes only; S5/S6 genes outside S2 excluded and counted',
        'human':'up and down separately, padj<0.05; background = nonmissing DESeq2 padj in each contrast',
        'singlecell':'top 50 cluster markers, log2FC>0.5, detected in >=25% of cluster; background = retained expressed genes',
        'interpretation':'Association with a gene set does not establish pathway activation, inhibition, direct binding, or causality.'},indent=2)+'\n')
    print(pd.DataFrame(audits).to_string(index=False))


if __name__ == '__main__':
    main()
