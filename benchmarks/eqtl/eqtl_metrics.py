#!/usr/bin/env python3
"""Report tissue-weighted eQTL sign AUROC and effect-size Spearman correlation."""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score


def weighted_metric(predictions, metric):
    scores, counts = [], []
    for _, tissue in predictions.groupby('tissue'):
        tissue = tissue.dropna(subset=['target', 'prediction'])
        if len(tissue) < 2:
            continue
        try:
            score = metric(tissue.target, tissue.prediction)
        except ValueError:
            continue
        if np.isfinite(score):
            scores.append(score)
            counts.append(len(tissue))
    return float(np.average(scores, weights=counts)) if counts else float('nan')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sign', help='Sign predictions Feather file')
    parser.add_argument('coefficient', help='Effect-size predictions Feather file')
    args = parser.parse_args()
    sign = pd.read_feather(args.sign)
    if 'metric_calculator' in sign:
        sign = sign[sign.metric_calculator.eq('TissueAggregatedAnnDataAUROC')]
    coefficient = pd.read_feather(args.coefficient)
    auroc = weighted_metric(sign, roc_auc_score)
    spearman = weighted_metric(coefficient, lambda truth, pred: spearmanr(truth, pred).statistic)
    print(f'sign tissue_weighted_mean_auroc: {auroc:.4f}')
    print(f'coefficient tissue_weighted_mean_spearmanr: {spearman:.4f}')


if __name__ == '__main__':
    main()
