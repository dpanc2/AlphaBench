# Enhancer–promoter interactions

Score how shuffling a 2 kb enhancer changes predicted K562 RNA-seq for its
target gene.

## Data

[`EP.GRCh38-c.tsv`](data/EP.GRCh38-c.tsv) is the prepared GRCh38 table for the
K562 CRISPRi enhancer–gene linking benchmark. It contains 10,356 candidate pairs,
including 471 positives. The runner uses `Regulated` as the label and the
`center` field as the zero-based start of the 2 kb replacement. The default
1 Mb window retains 6,606 pairs. See the
[ENCODE–rE2G study](https://doi.org/10.1038/s41586-026-10781-4) and the
[AlphaGenome enhancer–gene evaluation](https://doi.org/10.1038/s41586-025-10014-0).

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
evaluated when returned by the API.

## API rate limits

`--threads` controls concurrent requests and defaults to 5. If the API returns
a quota or rate-limit error, wait briefly, then rerun with `--threads 1` and
the same output directory. The runner reuses completed pair files and resumes
scoring the remaining pairs.
