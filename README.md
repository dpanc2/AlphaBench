# AlphaBench

AlphaBench runs independent AlphaGenome benchmarks using public experimental data.

```bash
git clone https://github.com/dpanc2/alphaBench.git AlphaBench
cd AlphaBench
```

| Benchmark | Start here |
| --- | --- |
| Gene expression | [Prepare ENCODE TPMs, predict expression, and compare](benchmarks/expression/README.md) |
| CAGI5 MPRA | [Run zero-shot DNase variant scoring](benchmarks/mpra/README.md) |
| ATAC-seq | [Prepare held-out regions and compare with ChromBPNet](benchmarks/atac/README.md) |
| eQTL | [Rebuild GTEx datasets, predict, and evaluate](benchmarks/eqtl/README.md) |
| Enhancer–promoter | [Score enhancer shuffles and calculate auPRC](benchmarks/pep/README.md) |

Run benchmark commands from the repository root. Downloaded inputs go under
`data/` and benchmark outputs under `results/`.

## Get an AlphaGenome API key

Follow the [AlphaGenome setup guide](https://www.alphagenomedocs.com/installation.html)
to create a free API key. In the terminal where you’ll run the benchmark, set
the environment variable used by AlphaBench:

```bash
export ALPHAGENOME_API_KEY='YOUR_API_KEY'
```

The hosted API runs inference remotely. Local AlphaGenome checkpoint inference
requires a CUDA GPU. Check the benchmark guide for its setup.

## Reference genomes

| Species | Assembly | Where it’s used |
| --- | --- | --- |
| Human | GRCh38 (hg38) | ATAC-seq, MPRA, and enhancer–promoter share `data/genomes/hg38.fa`. Human expression and eQTL use hg38 coordinates with the API. |
| Mouse | GRCm38 (mm10) | Mouse expression uses mm10 gene coordinates with the API and needs no local FASTA. |

Use UCSC chromosome names such as `chr1`, `chr2`, and `chrX`. See
[how to download the shared hg38 FASTA](data/README.md#shared-human-reference).

## Run with an agent

For agent-assisted runs, use the [benchmark skill](skills/benchmark/SKILL.md).

[License](LICENSE)
