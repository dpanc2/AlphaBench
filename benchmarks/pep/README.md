# Enhancer–promoter interactions

Score the included CRISPR enhancer–gene pairs using AlphaGenome ALL_FOLDS.
Each pair gets ten seeded 2 kb enhancer shuffles in a promoter-centered
1 Mb input. The score is the mean absolute effect on K562 total RNA-seq.

Use Python 3.11, CUDA-enabled JAX, access to the ALL_FOLDS weights, and the
[shared hg38 FASTA](../../data/README.md#shared-human-reference).

```bash
python -m pip install alphagenome==0.7.0 alphagenome-research==0.3.0 pysam scikit-learn
CUDA_VISIBLE_DEVICES=0 XLA_PYTHON_CLIENT_PREALLOCATE=false \
  python benchmarks/pep/run.py --position-offset 1
```

The runner evaluates the 6,606 pairs whose replacements fit within the
input. It uses the `Regulated` label and excludes missing scores from auPRC.
Find predictions, precision–recall curves, and `metrics.json` under
`results/pep/all_folds/`. Repeat the command to reuse completed pairs, or use
`--limit 2 --out results/pep/smoke` for a quick test.
