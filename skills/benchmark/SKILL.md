---
name: benchmark
description: Run or evaluate a model on alphaBench and report reproducible scores. Use for benchmark requests that name a model, such as Borzoi expression prediction.
---

# Run an alphaBench benchmark

1. From the repository root, read the benchmark guide linked in the root README. For Borzoi, use [the expression guide](../../benchmarks/borzoi/README.md) and the original held-out human test split.
2. For a demonstration, run one replicate on a bounded subset first. Use the full split when requested. Record the checkpoint, split, gene limit, target normalization, and whether inference was fresh or reused.
3. Check that predictions and truth share gene IDs and sample or ontology–assay IDs. For Borzoi, the runner exports sample accessions. The comparison script averages those into ontology–assay groups. Compare models on the same scored gene–group pairs.
4. Read the resulting `summary.json`. Report gene and cell Pearson correlations, matched gene and group counts, and the output directory. Label subset results as pilot scores. If inference or evaluation fails, report the error and the last completed stage.
