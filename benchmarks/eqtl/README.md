# eQTL

Rebuild GTEx SNV and indel datasets, score variants with AlphaGenome, and
evaluate effect size (Spearman) and direction (AUROC). Both metrics are
computed per tissue and weighted by the number of evaluated rows.

From the repository root:

```bash
bash scripts/setup.sh
source .venv/bin/activate
python -m pip install pyarrow scikit-learn
export ALPHAGENOME_API_KEY='YOUR_API_KEY'

python benchmarks/eqtl/prepare_datasets.py

dataset=indel_credible_sets_split
model=ALL_FOLDS
out="results/eqtl/fresh/$model/$dataset"
python benchmarks/eqtl/run_predictions.py \
  --variants "data/eqtl/rebuilt/$dataset.tsv" \
  --mapping data/eqtl/rebuilt/tissue_mapping.tsv \
  --model "$model" --out "$out"
python benchmarks/eqtl/eqtl_metrics.py \
  "$out/${dataset}_${model}_sign.feather" \
  "$out/${dataset}_${model}_coeff.feather"
```

Other datasets are `snv_eqtl` and `indel_credible_sets_merged`. Set
`model=FOLD_1` to evaluate that checkpoint.

SNV targets come from the [AlphaGenome evaluation data](https://storage.googleapis.com/alphagenome/evals/eqtl_variant_borzoi_coefficient_human_predictions.feather).
Indels come from the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/)
and retain credible-set variants with PIP ≥ 0.9.
