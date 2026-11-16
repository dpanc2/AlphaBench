#!/usr/bin/env python3
"""Run a strand- and assay-matched AlphaGenome 2 TPM benchmark."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import time
import json
import os
from pathlib import Path

import grpc
import numpy as np
import pandas as pd

from predict_expression import WINDOWS, read_genes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genes', required=True)
    parser.add_argument('--mapping', required=True, help='Ontology-assay to sample IDs JSON')
    parser.add_argument('--organism', choices=['human', 'mouse'], default='human')
    parser.add_argument('--window', choices=WINDOWS, default='128kb')
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--max-genes', type=int, default=32, help='0 uses all genes')
    parser.add_argument('--max-conditions', type=int, default=8, help='0 uses all matched conditions')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    if args.max_genes < 0 or args.max_conditions < 0:
        parser.error('Limits must be nonnegative')
    if args.workers < 1:
        parser.error('--workers must be positive')
    key = os.environ.get('ALPHAGENOME_API_KEY')
    if not key:
        parser.error('Set ALPHAGENOME_API_KEY')

    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, interval_scorers

    model = dna_client.create(key, timeout=30)
    organism = (dna_client.Organism.HOMO_SAPIENS if args.organism == 'human'
                else dna_client.Organism.MUS_MUSCULUS)
    width, mask_width = WINDOWS[args.window]
    scorer = interval_scorers.GeneMaskScorer(
        requested_output=dna_client.OutputType.RNA_SEQ,
        width=mask_width, aggregation_type=interval_scorers.IntervalAggregationType.MEAN,
    )
    genes = read_genes(args.genes, args.organism)
    if args.max_genes and len(genes) > args.max_genes:
        genes = genes.sample(args.max_genes, random_state=args.seed).sort_index()
    mapping = json.loads(Path(args.mapping).read_text())
    supported = model.supported_conditions(organism)
    observed = model.observed_metadata(organism).rna_seq
    matched = observed[
        observed.biosample_curie.isin(supported.loc[supported.rna_seq.eq('observed'), 'condition'])
        & observed.strand.isin(['+', '-'])
    ].copy()
    matched['group'] = matched.biosample_curie + ' ' + matched['Assay title']
    matched = matched[matched.group.isin(mapping)]
    conditions = sorted(matched.biosample_curie.unique())
    if not conditions:
        raise RuntimeError('No observed strand-specific AlphaGenome 2 tracks match the TPM assays')
    total_conditions = len(conditions)
    if args.max_conditions and len(conditions) > args.max_conditions:
        conditions = sorted(np.random.default_rng(args.seed).choice(
            conditions, args.max_conditions, replace=False
        ).tolist())
    matched = matched[matched.biosample_curie.isin(conditions)]
    args.out.mkdir(parents=True, exist_ok=True)
    plan = {
        'model': 'AlphaGenome 2', 'organism': args.organism, 'window': args.window,
        'genes': len(genes), 'conditions': conditions, 'matched_conditions_total': total_conditions,
        'assay_groups': int(matched.group.nunique()), 'seed': args.seed,
        'requests': int(sum(genes.gene_strand.eq(strand).sum() * count
                            for strand, count in matched.groupby('strand').biosample_curie.nunique().items())),
        'input_sha256': [hashlib.sha256(Path(path).read_bytes()).hexdigest()
                         for path in (args.genes, args.mapping)],
        'tracks_sha256': hashlib.sha256(matched.to_csv(index=False).encode()).hexdigest(),
    }
    plan_path = args.out / 'plan.json'
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        parser.error('Inputs or track selection changed. Use a new output directory.')
    plan_path.write_text(json.dumps(plan, indent=2) + '\n')
    matched.to_csv(args.out / 'matched_tracks.csv', index=False)
    print(json.dumps(plan), flush=True)
    if args.plan_only:
        return

    def score_gene(gene):
        interval = genome.Interval(gene.chrom, int(gene.TSS), int(gene.TSS),
                                   strand=gene.gene_strand).resize(width)
        row = {'gene_id': gene.gene_id}
        for condition, tracks in by_strand[gene.gene_strand]:
            for attempt in range(8):
                try:
                    scores = model.score_interval(interval, condition=condition,
                                                  organism=organism, interval_scorers=[scorer])[0]
                    break
                except grpc.RpcError as error:
                    if attempt == 7 or error.code() not in transient:
                        raise RuntimeError(
                            f'{gene.gene_id} / {condition}: {error.code().name}: {error.details()}'
                        ) from error
                    print(f'{error.code().name} for {gene.gene_id} / {condition}; retry {attempt + 1}', flush=True)
                    time.sleep(min(300, 10 * 2**attempt))
            if scores.n_obs == 0 or scores.n_vars == 0:
                continue
            indices = np.flatnonzero(scores.obs.gene_id.str.split('.').str[0]
                                     .eq(gene.gene_id_unversioned))
            if len(indices) != 1:
                continue
            for group, names in tracks.items():
                mask = scores.var.name.isin(names) & scores.var.strand.eq(gene.gene_strand)
                if mask.any():
                    values = np.asarray(scores.X[indices[0], mask.to_numpy()]).ravel()
                    finite = values[np.isfinite(values)]
                    if finite.size:
                        row[group] = float(finite.mean())
        return row

    transient = {grpc.StatusCode.UNAVAILABLE, grpc.StatusCode.DEADLINE_EXCEEDED,
                 grpc.StatusCode.RESOURCE_EXHAUSTED, grpc.StatusCode.INTERNAL}
    by_strand = {
        strand: [(condition, block.groupby('group')['name'].agg(list).to_dict())
                 for condition, block in matched[matched.strand.eq(strand)].groupby('biosample_curie')]
        for strand in ('+', '-')
    }
    checkpoint = args.out / 'scores.jsonl'
    records = {}
    if checkpoint.exists():
        with checkpoint.open('rb+') as handle:
            while line := handle.readline():
                if not line.endswith(b'\n'):
                    handle.truncate(handle.tell() - len(line))
                    break
                row = json.loads(line)
                records[row['gene_id']] = row
    remaining = [g for g in genes.itertuples() if g.gene_id not in records]
    print(f'Resuming: {len(records)}/{len(genes)} genes saved', flush=True)
    with checkpoint.open('a') as output, ThreadPoolExecutor(max_workers=args.workers) as pool:
        # Submit only one group of genes at a time so an API failure stops new work.
        for start in range(0, len(remaining), args.workers):
            for row in pool.map(score_gene, remaining[start:start + args.workers]):
                output.write(json.dumps(row) + '\n')
                output.flush()
                records[row['gene_id']] = row
                print(f"{len(records)}/{len(genes)} genes, {len(row) - 1} scores for {row['gene_id']}", flush=True)
    predictions = pd.DataFrame(list(records.values()))
    predictions.to_csv(args.out / 'predictions.csv', index=False)
    if len(predictions.columns) == 1:
        raise RuntimeError('No scores matched the observed assay and strand metadata')


if __name__ == '__main__':
    main()
