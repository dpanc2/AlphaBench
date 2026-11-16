"""Score ten enhancer shuffles per CRISPR enhancer–gene pair."""
import argparse
import dataclasses
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

K562 = 'EFO:0002067'
ASSAY = 'total RNA-seq'


def prepare_pairs(path, width, centered):
    data = pd.read_csv(path, sep='\t').drop_duplicates(
        ['measuredGeneEnsemblId', 'startTSS', 'chromStart'], keep='last'
    )
    enhancer = data.center + 1000
    centers = {
        'tss': data.startTSS,
        'mut': enhancer,
        'mid': (data.startTSS + enhancer) // 2,
    }
    starts = centers[centered] - width // 2
    ends = starts + width
    fits = ((data.chrom == data.chrTSS) & (starts < data.center)
            & (ends > data.center + 2000) & (starts <= data.startTSS)
            & (data.startTSS < ends))
    eligible = data.loc[fits].assign(interval_start=starts[fits], interval_end=ends[fits])
    return eligible.sort_values(['measuredGeneEnsemblId', 'startTSS', 'chromStart']), len(data)


def load_model(args):
    if args.version == 2:
        from alphagenome.models.v2 import dna_client, variant_scorers
        return dna_client.create(os.environ['ALPHAGENOME_API_KEY'], timeout=30), variant_scorers
    if args.backend == 'local':
        from alphagenome.models import variant_scorers
        from alphagenome_research.model import dna_model

        settings = dict(dna_model.default_organism_settings())
        human = dna_model.Organism.HOMO_SAPIENS
        settings[human] = dataclasses.replace(settings[human], fasta_path=args.genome,
                                             calibration_path=None, pas_feather_path=None)
        settings[dna_model.Organism.MUS_MUSCULUS] = dna_model.OrganismSettings()
        return dna_model.create_from_huggingface('all_folds', organism_settings=settings), variant_scorers
    try:
        from alphagenome.models.v1 import dna_client, variant_scorers
    except ImportError:
        from alphagenome.models import dna_client, variant_scorers
    return dna_client.create(os.environ['ALPHAGENOME_API_KEY'],
                             model_version=dna_client.ModelVersion.ALL_FOLDS,
                             timeout=30), variant_scorers


def score_pair(pair, index, fasta, model, scorers, args):
    from alphagenome.data import genome

    reference = fasta.fetch(pair.chrom, pair.center, pair.center + 2000).upper()
    if len(reference) != 2000:
        raise ValueError(f'Incomplete hg38 sequence for {pair.chrom}:{pair.center}')
    rng = np.random.default_rng(np.random.SeedSequence([args.seed, index]))
    interval = genome.Interval(pair.chrTSS, pair.interval_start, pair.interval_end)
    variants = [
        genome.Variant(chromosome=pair.chrom, position=pair.center + args.position_offset,
                       reference_bases=reference,
                       alternate_bases=''.join(rng.permutation(list(reference))))
        for _ in range(10)
    ]
    options = {'conditions': [K562] * 10} if args.version == 2 else {}
    responses = model.score_variants(
        [interval] * 10, variants, variant_scorers=[scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
        max_workers=args.threads, progress_bar=False, **options,
    )
    rows = []
    for shuffle, response in enumerate(responses):
        scores = scorers.tidy_scores([response])
        if scores is None or scores.empty:
            continue
        ontology = 'biosample_curie' if 'biosample_curie' in scores else 'ontology_curie'
        matches = scores.loc[
            scores.gene_id.str.split('.').str[0].eq(pair.measuredGeneEnsemblId)
            & scores[ontology].eq(K562) & scores['Assay title'].eq(ASSAY)
        ]
        if matches.empty:
            continue
        row = {'shuffle': shuffle, 'raw_score': matches.raw_score.mean()}
        if 'quantile_score' in matches:
            row['quantile_score'] = matches.quantile_score.mean()
        rows.append(row)
    return pd.DataFrame(rows, columns=['shuffle', 'raw_score', 'quantile_score'])


def mean_abs(scores, column):
    if len(scores) != 10 or column not in scores or not np.isfinite(scores[column]).all():
        return np.nan
    return float(scores[column].abs().mean())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome', default='data/genomes/hg38.fa')
    parser.add_argument('--pairs', default=Path(__file__).parent / 'data/EP.GRCh38-c.tsv', type=Path)
    parser.add_argument('--out', default=Path('results/pep/api_v1'), type=Path)
    parser.add_argument('--version', choices=[1, 2], type=int, default=1)
    parser.add_argument('--backend', choices=['api', 'local'], default='api')
    parser.add_argument('--centered', choices=['tss', 'mut', 'mid'], default='tss')
    parser.add_argument('--window', type=int, choices=[16384, 131072, 524288, 1048576], default=1048576)
    parser.add_argument('--threads', type=int, default=5)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--position-offset', type=int, choices=[0, 1], default=1)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if args.threads < 1 or (args.limit is not None and args.limit < 1):
        parser.error('--threads and --limit must be positive')
    if args.version == 2 and args.backend == 'local':
        parser.error('AlphaGenome 2 requires the API backend')
    eligible, total = prepare_pairs(args.pairs, args.window, args.centered)
    print(f'{len(eligible)}/{total} pairs fit the selected window', flush=True)
    if eligible.empty:
        parser.error('No pairs fit the selected window')
    if args.prepare_only:
        args.out.mkdir(parents=True, exist_ok=True)
        eligible.to_csv(args.out / 'eligible_pairs.tsv', sep='\t', index=False)
        return
    if args.backend == 'api' and not os.environ.get('ALPHAGENOME_API_KEY'):
        parser.error('Set ALPHAGENOME_API_KEY before scoring')

    import pysam

    fasta_path = Path(args.genome).resolve(strict=True)
    stat = fasta_path.stat()
    settings = {'backend': args.backend, 'version': args.version, 'centered': args.centered,
                'window': args.window, 'seed': args.seed, 'position_offset': args.position_offset,
                'pairs_sha256': hashlib.sha256(args.pairs.read_bytes()).hexdigest(),
                'genome': str(fasta_path), 'genome_size': stat.st_size, 'genome_mtime_ns': stat.st_mtime_ns}
    args.out.mkdir(parents=True, exist_ok=True)
    identity = args.out / 'settings.json'
    if identity.exists() and json.loads(identity.read_text()) != settings:
        parser.error('Run settings changed. Choose a new --out directory.')
    cache = args.out / 'pairs'
    cache.mkdir(exist_ok=True)
    identity.write_text(json.dumps(settings, indent=2) + '\n')
    eligible.to_csv(args.out / 'eligible_pairs.tsv', sep='\t', index=False)

    if not Path(args.genome + '.fai').exists():
        pysam.faidx(args.genome)
    predictions = []
    model = scorers = None
    selected = eligible.head(args.limit) if args.limit else eligible
    with pysam.FastaFile(args.genome) as fasta:
        missing = set(selected.chrom) - set(fasta.references)
        if missing:
            parser.error(f'FASTA is missing chromosomes: {", ".join(sorted(missing))}')
        for index, pair in enumerate(selected.itertuples()):
            path = cache / f'{index:05d}.tsv'
            if path.exists():
                scores = pd.read_csv(path, sep='\t')
            else:
                if model is None:
                    model, scorers = load_model(args)
                scores = score_pair(pair, index, fasta, model, scorers, args)
                temporary = path.with_suffix('.tmp')
                scores.to_csv(temporary, sep='\t', index=False)
                temporary.replace(path)
                print(f'{index + 1}/{len(selected)} {pair.measuredGeneEnsemblId}: {len(scores)}/10 scores',
                      flush=True)
            predictions.append({
                'gene_id': pair.measuredGeneEnsemblId, 'chrom': pair.chrom, 'center': pair.center,
                'distance': pair.distanceToTSS, 'label': str(pair.Regulated).upper() == 'TRUE',
                'score': mean_abs(scores, 'raw_score'),
                'quantile_score': mean_abs(scores, 'quantile_score'),
            })
    pd.DataFrame(predictions).to_csv(args.out / 'predictions.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
