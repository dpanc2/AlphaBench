# Borzoi expression demo

Run human Borzoi replicate 0 on the first 32 rows of the held-out test split, then
compare the predicted RNA-seq signal with measured TPMs. The runner adapts
the [GENA_LM Borzoi workflow](https://github.com/dpanc2/GENA_LM/tree/26d81ed41d0fb24f0ef97b92b414778c8f8816d8/downstream_tasks/expression_prediction/benchmarks/borzoi).

Use Linux with a CUDA GPU, Python 3.11, CUDA-enabled PyTorch, and
`borzoi-pytorch==0.4.4`. Install the [shared hg38 FASTA](../mpra/README.md#shared-genome).
Obtain the original TPM matrix or prepare one with the
[expression guide](../expression/README.md#prepare-your-own-ground-truth).
Use a separate output directory for each run.

From the repository root:

```bash
mkdir -p data/borzoi
source_url=https://raw.githubusercontent.com/dpanc2/GENA_LM/26d81ed41d0fb24f0ef97b92b414778c8f8816d8/downstream_tasks/expression_prediction/benchmarks/borzoi
for file in targets_human.txt qnorm_id_to_targets_map.csv \
    gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz \
    human.test.forward.csv human.test.reverse.csv; do
  curl --fail --location "$source_url/$file" -o "data/borzoi/$file" || break
done

python benchmarks/borzoi/predict.py borzoi test --replicates 0 --max-rows 32 \
  --borzoi-dir data/borzoi --device cuda:0 --output-dir results/borzoi_demo
python scripts/compare_expression.py \
  --predictions results/borzoi_demo/borzoi_predictions_test_replicate0_pytorch.csv \
  --truth data/human_id_full_tpms.csv \
  --mapping benchmarks/expression/data/alphagenome_track_to_qnorm_ids_all.json \
  --out results/borzoi_demo/raw_tpm
```

The model uses 524,288 bp input and gene-strand-matched exon overlap.
The published checkpoint conversion is `johahi/borzoi-replicate-0`.
Omit `--max-rows` for the full split. The evaluation writes
`raw_tpm/summary.json` and per-gene and per-group correlations.
Use a qnorm TPM matrix as `--truth` for the normalized comparison.
