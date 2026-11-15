# CAGI5 MPRA

Score CAGI5 variants with AlphaGenome ALL_FOLDS and compare predicted DNase
effects with reporter activity. Enformer and Borzoi specify the tissue
mappings. Both use a 1 Mb input and the same 501 bp ALT-minus-REF scorer.

Use Python 3.11+, CUDA-enabled JAX, access to the ALL_FOLDS weights, and the
[shared hg38 FASTA](../../data/README.md#shared-human-reference). From the root:

```bash
python -m pip install alphagenome==0.7.0 alphagenome-research==0.3.0
bash benchmarks/mpra/setup.sh
mkdir -p data/mpra
curl --fail --location 'https://zenodo.org/records/17453603/files/Kircher_GRCh38_ALL.tsv?download=1' \
  -o data/mpra/Kircher_GRCh38_ALL.tsv
python benchmarks/mpra/prepare.py
config=results/mpra/cagi5.filtered.yaml
alphabench download --config "$config"
alphabench process --config "$config"
CUDA_VISIBLE_DEVICES=0 alphabench run --config "$config"
```

Preparation removes one observation whose reference allele differs from hg38.
The run evaluates 29,034 observations across 22 assay groups per mapping.
Find predictions and Pearson and Spearman correlations under `results/mpra/`.
To resume, add `--resume <run-directory>` to the run command.
