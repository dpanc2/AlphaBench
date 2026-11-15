#!/usr/bin/env python3
"""Predict GTEx eQTL effects with resumable AlphaGenome variant scoring."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import time

import grpc
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from tqdm import tqdm

DATASETS = ('snv_eqtl', 'indel_credible_sets_split', 'indel_credible_sets_merged')
# Preserve caches made by the previous runner with the same scoring protocol.
SCORING_ID = '2c7f009fb09ef3130a0c6ec538738c8994ebde83dbc076233b959d2d9558d92f'
RETRYABLE = {grpc.StatusCode.UNAVAILABLE, grpc.StatusCode.DEADLINE_EXCEEDED,
             grpc.StatusCode.RESOURCE_EXHAUSTED}


def score_variant(model, variant_id):
    chrom, position, ref, alt = variant_id.rsplit('_', 3)
    variant = genome.Variant(chrom, int(position), ref, alt)
    interval = variant.reference_interval.resize(
        dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']
    )
    interval = genome.Interval(interval.chromosome, interval.start, interval.end, strand='+')
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']
    for attempt in range(10):
        try:
            score = model.score_variant(interval, variant, [scorer])[0]
            break
        except grpc.RpcError as error:
            if attempt == 9 or error.code() not in RETRYABLE:
                raise
            time.sleep(60)
    return variant_scorers.tidy_scores(
        [[score]], match_gene_strand=True, include_extended_metadata=True
    )


def export_predictions(rows, variants, destination):
    keys = ['variant_id', 'gene_id', 'tissue']
    targets = variants.assign(
        variant_id=variants.variant + '_b38', tissue=variants.alphagenome_label
    )[keys + ['beta']]
    table = pd.DataFrame(rows, columns=keys + ['prediction'])
    table['variant_id'] += '_b38'
    table = table.merge(targets, on=keys, how='left').sort_values('tissue', kind='stable')
    table['prediction'] = table.prediction.astype('float32')
    table['variant_scorer'] = 'gene_exons_contained_diff_mean;OutputType.RNA_SEQ;log_0001p_after_agg'
    table['output_type'] = 'RNA_SEQ'
    table['metric_name'] = None
    columns = keys + ['prediction', 'target', 'variant_scorer', 'output_type',
                      'metric_calculator', 'metric_name']
    coefficient = table.assign(target=table.beta.astype('float32'),
                               metric_calculator='TissueAggregatedAnnDataCorrelations')
    labels = (table.beta >= 0).astype('float32').where(table.beta.notna())
    sign = pd.concat([table.assign(target=labels, metric_calculator=metric)
                      for metric in ('TissueAggregatedAnnDataAUROC', 'TissueAggregatedAnnDataAUPRC')])
    for suffix, frame in [('coeff', coefficient), ('sign', sign)]:
        frame[columns].reset_index(drop=True).to_feather(f'{destination}_{suffix}.feather')


def run(variants, mapping, model, out, cache, limit):
    variants, mapping, out, cache = map(Path, (variants, mapping, out, cache))
    inputs = pd.read_csv(variants, sep='\t')
    if limit:
        inputs = inputs[inputs.variant.isin(inputs.variant.drop_duplicates().head(limit))]
    if inputs.empty:
        raise ValueError('No variants selected')
    tissues = pd.read_csv(mapping, sep='\t')[['dataset_id', 'alphagenome_label']]
    pairs = inputs.merge(tissues, on='dataset_id', how='left', validate='many_to_one')
    if pairs.alphagenome_label.isna().any():
        raise ValueError('Input contains dataset IDs absent from the tissue mapping')

    sdk = importlib.metadata.version('alphagenome')
    serialized = inputs.to_csv(sep='\t', index=False).encode()
    manifest = {'model': model, 'input_sha256': hashlib.sha256(serialized).hexdigest(),
                'mapping_sha256': hashlib.sha256(mapping.read_bytes()).hexdigest(),
                'source_sha256': SCORING_ID, 'sdk': sdk}
    out.mkdir(parents=True, exist_ok=True)
    cache = cache / model
    cache.mkdir(parents=True, exist_ok=True)
    manifest_path = out / 'inputs.json'
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        if any(previous.get(key) != value for key, value in manifest.items()):
            raise ValueError('Input or settings changed. Choose a new output folder.')
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    (out / variants.name).write_bytes(serialized)

    all_genes = {}
    # Reuse the same cache across the standard datasets, including both indel views.
    for path in {variants, *(variants.parent / f'{name}.tsv' for name in DATASETS)}:
        if path.is_file():
            table = pd.read_csv(path, sep='\t', usecols=['variant', 'gene_id'])
            for variant_id, genes in table.groupby('variant').gene_id:
                all_genes.setdefault(variant_id, set()).update(genes)

    client = None
    rows = []
    groups = pairs.groupby('variant', sort=False).gene_id
    for variant_id, genes in tqdm(groups, total=groups.ngroups):
        gene_ids = sorted(all_genes[variant_id])
        identity = json.dumps([variant_id, gene_ids, sdk, SCORING_ID], sort_keys=True)
        path = cache / (hashlib.sha256(identity.encode()).hexdigest() + '.parquet')
        if path.is_file():
            tidy = pd.read_parquet(path)
        else:
            if client is None:
                client = dna_client.create(os.environ['ALPHAGENOME_API_KEY'],
                                           model_version=dna_client.ModelVersion[model])
            tidy = score_variant(client, variant_id)
            tidy = tidy.loc[tidy.gene_id.str.split('.').str[0].isin(gene_ids),
                            ['gene_id', 'gtex_tissue', 'raw_score']].copy()
            temporary = path.with_suffix('.tmp')
            tidy.to_parquet(temporary, index=False)
            temporary.replace(path)
        scores = tidy.groupby([tidy.gene_id.str.split('.').str[0], 'gtex_tissue']).raw_score.mean()
        rows.extend((variant_id, gene, tissue, scores.get((gene, tissue), 0.0))
                    for gene in genes.unique() for tissue in tissues.alphagenome_label)

    export_predictions(rows, pairs, out / f'{variants.stem}_{model}')
    (out / 'complete.json').write_text(json.dumps({
        'status': 'complete', 'model': model, 'variants': groups.ngroups, 'pairs': len(pairs)
    }, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variants', required=True)
    parser.add_argument('--mapping', required=True)
    parser.add_argument('--model', choices=['ALL_FOLDS', 'FOLD_1'], required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--cache', default='results/eqtl/fresh/cache')
    parser.add_argument('--limit', type=int, help='Only score the first N unique variants')
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error('--limit must be positive')
    run(**vars(args))
