"""Validate public inputs; preserve species and original study comparisons."""
from pathlib import Path
import gzip
import hashlib
import json
import re
import tarfile
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw'
OUT = ROOT / 'data/processed'
OUT.mkdir(parents=True, exist_ok=True)


def samples(path):
    text = gzip.open(path, 'rt').read()
    records = []
    for block in text.split('^SAMPLE = ')[1:]:
        record = {'accession': block.splitlines()[0]}
        for line in block.splitlines()[1:]:
            if line.startswith('!Sample_') and ' = ' in line:
                key, value = line.split(' = ', 1)
                record.setdefault(key.removeprefix('!Sample_'), []).append(value)
        records.append(record)
    return records


def main():
    book = RAW / 'PhysiolGen_SupplementaryTables_2019.xls'
    metadata = json.loads((RAW / 'figshare_8855588.json').read_text())
    info = next(x for x in metadata['files'] if x['name'] == book.name)
    assert hashlib.md5(book.read_bytes()).hexdigest() == info['computed_md5']
    tables = {}
    for number, expected in [(2, 617), (5, 258), (6, 107)]:
        t = pd.read_excel(book, sheet_name=f'Supplementary TABLE S{number}', skiprows=1)
        assert t.shape == (expected, 4)
        t.columns = ['symbol', 'ensembl_id', 'description', 'reported_signed_fold_change']
        assert t.ensembl_id.is_unique and t.ensembl_id.str.startswith('ENSRNOG').all()
        assert np.isfinite(t.reported_signed_fold_change).all()
        # Preserve the source value. In particular, S5 Ccl20/Mmp9 signs differ
        # from the narrative description; do not silently reverse them.
        t['species'] = 'Rattus norvegicus'
        t['source_table'] = f'S{number}'
        t['comparison_as_reported'] = {2:'PIA vs naive',5:'RTD-1 vs saline day 3',6:'RTD-1 vs saline day 5'}[number]
        t.to_csv(OUT / f'rat_2019_S{number}.tsv', sep='\t', index=False)
        tables[number] = t
    union = set(tables[5].ensembl_id) | set(tables[6].ensembl_id)
    summary = {'arthritis_genes_S2':617, 'day3_genes_S5':258, 'day5_genes_S6':107,
               'day3_day5_union':len(union), 'day3_day5_overlap':len(set(tables[5].ensembl_id)&set(tables[6].ensembl_id)),
               'day1_gene_table_available':False,
               'day3_genes_outside_S2': sorted(set(tables[5].ensembl_id)-set(tables[2].ensembl_id)),
               'day5_genes_outside_S2': sorted(set(tables[6].ensembl_id)-set(tables[2].ensembl_id)),
               'direction_caution':'S5 Ccl20 and Mmp9 have positive fold changes despite narrative suppression. Values retained verbatim; no reversal claim.'}
    (OUT / 'rat_2019_audit.json').write_text(json.dumps(summary, indent=2)+'\n')
    counts = pd.read_csv(RAW / 'GSE278567_Partek_THP1_Gene_counts.txt.gz', sep='\t')
    assert counts.gene_id.is_unique
    rows = []
    for s in samples(RAW / 'GSE278567_family.soft.gz'):
        label = s['description'][0]
        title = s['title'][0]
        assert label in counts.columns
        dose_match = re.search(r'(1|3|10) [µμ]g/ml', title)
        dose = dose_match.group(1) if dose_match else '0'
        lps = 'LPS' in title
        rows.append({'sample':label, 'accession':s['accession'], 'title':title,
                     'species':s['organism_ch1'][0], 'dose_ug_ml':int(dose),
                     'lps':lps, 'group':('LPS_' if lps else 'basal_')+dose,
                     'replicate_label':re.search(r'exp(\d+)', title).group(1)})
    design = pd.DataFrame(rows)
    assert len(design) == 24 and design['sample'].is_unique
    assert (design.groupby('group').size() == 3).all()
    assert design.species.eq('Homo sapiens').all()
    x = counts[design['sample']].to_numpy(float)
    assert np.isfinite(x).all() and (x >= 0).all()
    # E/M allocation yields fractional estimated counts, not normalized values.
    # Explicit nearest-integer approximation for DESeq2; never use the TMM file.
    rounded = np.floor(x + 0.5).astype(np.int64)
    pd.DataFrame(rounded, index=counts.gene_id, columns=design['sample']).rename_axis('gene_id').to_csv(OUT / 'human_counts_rounded.tsv', sep='\t')
    counts[['gene_id','gene_name','gene_type']].to_csv(OUT / 'human_genes.tsv', sep='\t', index=False)
    design.to_csv(OUT / 'human_samples.tsv', sep='\t', index=False)
    qc = pd.DataFrame({'sample':design['sample'], 'estimated_library_size':x.sum(axis=0),
                       'rounded_library_size':rounded.sum(axis=0), 'detected_genes':(x>0).sum(axis=0),
                       'fraction_entries_noninteger':(np.abs(x-np.rint(x))>1e-6).mean(axis=0)})
    qc.to_csv(OUT / 'human_count_qc.tsv', sep='\t', index=False)
    sc_rows = []
    for s in samples(RAW / 'GSE269120_family.soft.gz'):
        sc_rows.append({'accession':s['accession'], 'title':s['title'][0],
                        'species':s['organism_ch1'][0], 'characteristics':'; '.join(s.get('characteristics_ch1',[]))})
    assert len(sc_rows) == 2
    pd.DataFrame(sc_rows).to_csv(OUT / 'rat_singlecell_samples.tsv', sep='\t', index=False)
    with tarfile.open(RAW / 'GSE269120_RAW.tar') as archive:
        for member in archive:
            assert member.isfile() and '/' not in member.name and member.name.endswith('.h5')
            target = RAW / member.name
            if not target.exists():
                target.write_bytes(archive.extractfile(member).read())
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
