# Data

The bundled `file_mappings_borzoi_and_v1_human.csv` and
`file_mappings_borzoi_all.csv` provide ENCODE sample IDs for the
[expression workflow](../benchmarks/expression/README.md).

## Shared human reference

ATAC-seq, MPRA, and enhancer–promoter use `data/genomes/hg38.fa`.
The ATAC runner downloads it if needed. To set it up separately:

```bash
mkdir -p data/genomes
curl --fail --location --retry 3 \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
  -o data/genomes/hg38.fa.gz
gzip -dc data/genomes/hg38.fa.gz > data/genomes/hg38.fa
```
