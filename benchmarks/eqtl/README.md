# AlphaGenome eQTL benchmark

## Layout

| File or directory | Contents |
| --- | --- |
| `benchmarks/eqtl/prepare_datasets.py` | Builds the three benchmark sets |
| `benchmarks/eqtl/data/tissue_mapping.tsv` | Maps dataset IDs to AlphaGenome tissues |
| `benchmarks/eqtl/run_predictions.py` | Generates AlphaGenome predictions |
| `benchmarks/eqtl/eqtl_metrics.py` | Calculates benchmark metrics |
| `data/eqtl/rebuilt/` | Generated SNV and indel sets |
| `results/eqtl/fresh/` | Generated predictions and scoring cache |

## Prepare datasets

Run from the repository root:

```bash
bash scripts/setup.sh
source .venv/bin/activate
python -m pip install pyarrow scikit-learn
python benchmarks/eqtl/prepare_datasets.py
```

Preparation writes `snv_eqtl.tsv`, `indel_credible_sets_merged.tsv`, and
`indel_credible_sets_split.tsv`. SNV targets come from the
[AlphaGenome evaluation data](https://storage.googleapis.com/alphagenome/evals/eqtl_variant_borzoi_coefficient_human_predictions.feather).
Indels come from the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/)
with PIP ≥ 0.9. The split set keeps the two cerebellum datasets separate.

## Predictions

```bash
export ALPHAGENOME_API_KEY='YOUR_API_KEY'
dataset=snv_eqtl
model=ALL_FOLDS
out="results/eqtl/fresh/$model/$dataset"
python benchmarks/eqtl/run_predictions.py \
  --variants "data/eqtl/rebuilt/$dataset.tsv" \
  --mapping data/eqtl/rebuilt/tissue_mapping.tsv \
  --model "$model" --out "$out"
```

Choose any prepared set. The required input columns are `gene_id`, `variant`,
`beta`, and `dataset_id`. The runner supports `ALL_FOLDS` and `FOLD_1`.
It writes `<set>_<model>_coeff.feather` and `<set>_<model>_sign.feather`
to the output directory.

## Metrics

```bash
python benchmarks/eqtl/eqtl_metrics.py \
  "$out/${dataset}_${model}_sign.feather" \
  "$out/${dataset}_${model}_coeff.feather"
```

The script prints tissue-weighted AUROC for effect direction and Spearman
correlation for effect size.
