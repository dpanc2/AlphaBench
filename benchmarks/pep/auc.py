"""Calculate auPRC from saved enhancer–promoter pair scores."""
import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.metrics import auc, precision_recall_curve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path('results/pep/api_v1'))
    args = parser.parse_args()
    predictions = pd.read_csv(args.out / 'predictions.tsv', sep='\t')
    bins = [('all', 0, float('inf')), ('0-6kb', 0, 6000), ('0-10kb', 0, 10000),
            ('10-100kb', 10000, 100000), ('100-2500kb', 100000, 2500000)]
    metrics = []
    for name, start, end in bins:
        candidates = predictions[predictions.distance.between(start, end, inclusive='left')]
        for column in ('score', 'quantile_score'):
            if column not in predictions or (column == 'quantile_score' and predictions[column].isna().all()):
                continue
            subset = candidates.dropna(subset=[column])
            value = None
            if subset.label.nunique() == 2:
                precision, recall, _ = precision_recall_curve(subset.label, subset[column])
                value = float(auc(recall, precision))
                suffix = name if column == 'score' else f'{name}_quantile'
                pd.DataFrame({'precision': precision, 'recall': recall}).to_csv(
                    args.out / f'pr_{suffix}.csv', index=False
                )
            metrics.append({'distance': name, 'score_type': column, 'n_pairs': len(subset),
                            'n_missing': len(candidates) - len(subset),
                            'n_positive': int(subset.label.sum()), 'auprc': value})
    (args.out / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
