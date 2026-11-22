"""Download public Objective 1 inputs and record SHA256 provenance.

Existing complete files are reused. Downloads use .part files and atomic rename;
passing --manifest-only records the already downloaded inputs without networking.
"""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import json
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw'
SOURCES = {
    'figshare_8855588.json': 'https://api.figshare.com/v2/articles/8855588',
    'PhysiolGen_SupplementaryTables_2019.xls':'https://ndownloader.figshare.com/files/16224056',
    'PhysiolGen_SupplementaryFigures_2019.pdf':'https://ndownloader.figshare.com/files/16224050',
    'paper_2019.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6962595/',
    'GSE278567_Partek_THP1_Gene_counts.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278567/suppl/GSE278567_Partek_THP1_Gene_counts.txt.gz',
    'GSE278567_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278567/soft/GSE278567_family.soft.gz',
    'GSE269120_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269120/soft/GSE269120_family.soft.gz',
    'GSE269120_RAW.tar':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269120/suppl/GSE269120_RAW.tar',
    'Ensembl2Reactome.txt':'https://reactome.org/download/current/Ensembl2Reactome.txt',
    'reactome_version.txt':'https://reactome.org/ContentService/data/database/version',
    'paper_2017.html':'https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0187868',
    'paper_mtd_2021.html':'https://www.nature.com/articles/s41598-021-02619-y',
    'MTD12813_2021_supplement.pdf':'https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-021-02619-y/MediaObjects/41598_2021_2619_MOESM1_ESM.pdf',
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest-only', action='store_true')
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, url in SOURCES.items():
        target = RAW / name
        if not target.exists():
            if args.manifest_only:
                raise FileNotFoundError(target)
            temporary = target.with_suffix(target.suffix+'.part')
            with urllib.request.urlopen(url, timeout=120) as response, temporary.open('wb') as sink:
                while chunk := response.read(1024*1024):
                    sink.write(chunk)
            temporary.replace(target)
        digest = hashlib.sha256()
        with target.open('rb') as stream:
            while chunk := stream.read(1024*1024): digest.update(chunk)
        rows.append({'file':'data/raw/'+name,'url':url,'bytes':target.stat().st_size,
                     'sha256':digest.hexdigest(),
                     'local_mtime_utc':dt.datetime.fromtimestamp(target.stat().st_mtime,dt.timezone.utc).isoformat()})
    (ROOT / 'data/source_manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(f'Validated presence and checksummed {len(rows)} source files')


if __name__ == '__main__':
    main()
