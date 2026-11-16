# Enhancer–promoter interactions

Score how shuffling a 2 kb enhancer changes predicted K562 RNA-seq for its
target gene. The CRISPR enhancer–gene pairs are included in
`data/EP.GRCh38-c.tsv`.

## Run

From the repository root, use Python 3.11+ and the
[shared hg38 FASTA](../../data/README.md#shared-human-reference) at
`data/genomes/hg38.fa`. Chromosome names should be `chr1`, `chr2`, …, `chrX`.

```bash
python -m pip install alphagenome==2.0.0 pysam scikit-learn
export ALPHAGENOME_API_KEY='YOUR_API_KEY'
python benchmarks/pep/run.py
python benchmarks/pep/auc.py
```

The default is AlphaGenome 1 ALL_FOLDS with 1 Mb windows centered at the TSS.
Use `--centered mid` to center between the TSS and enhancer or `--centered mut`
to center on the enhancer. Both the TSS and the full replacement must fit
within the window. To try AlphaGenome 2, use `--version 2 --out results/pep/api_v2`
and pass the same `--out` to `auc.py`.

Scores and auPRC are saved in `results/pep/api_v1/`. The metrics include
the subset under 6 kb and the broader distance groups. Quantile scores are
evaluated when returned by the API. Rerun the same command to reuse cached
pairs. If the API rate-limits the run, lower `--threads` and retry.
