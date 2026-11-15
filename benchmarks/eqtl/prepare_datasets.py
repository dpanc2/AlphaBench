#!/usr/bin/env python3
"""Rebuild SNV targets from official evaluation data and indels from eQTL Catalogue."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import pandas as pd
import requests


def prepare(out):
    mapping_path = Path(__file__).with_name("data") / "tissue_mapping.tsv"
    out = Path(out)
    raw = out / 'raw'; raw.mkdir(parents=True, exist_ok=True)
    urls = {}
    def fetch(url, name):
        path = raw / name; urls[name] = url
        if not path.exists():
            for attempt in range(3):
                try:
                    r = requests.get(url, timeout=180); r.raise_for_status()
                    temp = path.with_suffix('.tmp'); temp.write_bytes(r.content); temp.replace(path)
                    break
                except requests.RequestException:
                    if attempt == 2: raise
        return path
    mapping = pd.read_csv(mapping_path, sep='\t')
    snv_path = fetch('https://storage.googleapis.com/alphagenome/evals/eqtl_variant_borzoi_coefficient_human_predictions.feather', 'snv_coeff.feather')
    snv = pd.read_feather(snv_path).dropna(subset=['target']).copy()
    snv['gene_id'] = snv.gene_id.str.split('.').str[0]
    snv['variant'] = snv.variant_id.str.removesuffix('_b38')
    snv = snv.merge(mapping, left_on='tissue', right_on='alphagenome_label', validate='many_to_one')
    snv = snv.rename(columns={'target':'beta'})[['gene_id','variant','beta','dataset_id']]
    gtf = pd.read_feather(fetch('https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather', 'gencode_v46.feather'))
    genes = set(gtf.gene_id.dropna().str.split('.').str[0])
    def tissue(row):
        name = f'{row.dataset_id}.credible_sets.tsv.gz'
        path = fetch(f'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000015/{row.dataset_id}/{name}', name)
        f = pd.read_csv(path, sep='\t')
        alleles = f.variant.str.rsplit('_', n=3, expand=True)
        f = f[(f.pip >= .9) & ((alleles[2].str.len()!=1) | (alleles[3].str.len()!=1)) & f.gene_id.isin(genes)].copy()
        f = f.drop_duplicates(['variant','gene_id'])
        f['dataset_id'] = row.dataset_id; f['tissue'] = row.tissue_label
        print(row.dataset_id, len(f), flush=True)
        return f
    with ThreadPoolExecutor(max_workers=4) as pool:
        split = pd.concat(pool.map(tissue, mapping.itertuples()), ignore_index=True)
    # Keep the smallest p-value per variant, gene and tissue label.
    merged = split.sort_values('pvalue', kind='stable').drop_duplicates(['variant','gene_id','tissue'])
    tables = {'snv_eqtl':snv,'indel_credible_sets_split':split,'indel_credible_sets_merged':merged}
    audit = {}
    for name, table in tables.items():
        table.to_csv(out / f'{name}.tsv', sep='\t', index=False)
        audit[name] = {'pairs': len(table), 'unique_variants': table.variant.nunique()}
    shutil.copy2(mapping_path, out / 'tissue_mapping.tsv')
    (out/'preparation_manifest.json').write_text(json.dumps({'datasets':audit,'raw_files':{name:{'url':url,'sha256':hashlib.sha256((raw/name).read_bytes()).hexdigest()} for name,url in sorted(urls.items())}},indent=2)+'\n')
    print(json.dumps(audit,indent=2),flush=True)

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',default='data/eqtl/rebuilt')
    prepare(**vars(p.parse_args()))
