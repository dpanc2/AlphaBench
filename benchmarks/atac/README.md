# ATAC-seq benchmark

Compare AlphaGenome FOLD_0, FOLD_1 and ChromBPNet on IMR-90 ATAC-seq data.
Requires Linux, Python 3.10/3.11 and an AlphaGenome API key. The bundled ChromBPNet setup runs on CPU.

## 1. Open the project and load your key

From the repository root:

```bash
bash scripts/setup.sh
source .venv/bin/activate
export ALPHAGENOME_API_KEY='paste-your-key-here'
```

## 2. Prepare the BED files once

Open [`prepare_inputs.ipynb`](prepare_inputs.ipynb), select the `.venv`
kernel and run all cells. It removes peaks overlapping AlphaGenome training
regions expanded to 1 Mb and saves `fold0_test_regions.bed` (11,860 peaks)
and `fold1_test_regions.bed` (13,812 peaks) in `data/atac/`.
All three models run on the latter set. FOLD_0 is a diagnostic because
most of those regions overlap its training partitions.
If needed, install its dependencies first:

```bash
pip install alphagenome==0.7.0 pyranges==0.1.4 pandas
```

## 3. Try two regions first

```bash
bash benchmarks/atac/run_benchmark.sh results/atac/quick_test 2
```

This installs the benchmark dependencies, downloads missing signal/model files
and tests all three models on two regions. Prepare the BEDs in the notebook first.
A small results table appears when it finishes. These two-region metrics are
only a test, not the final benchmark results.

## 4. Run the full benchmark

```bash
mkdir -p results/atac
nohup bash benchmarks/atac/run_benchmark.sh results/atac/full_run \
  > results/atac/full_run.log 2>&1 &
```

Follow progress in `results/atac/full_run.log`. Each model also writes a log in its output directory.

## Results and restarting

The final table is `results/atac/full_run/summary.csv`. Each model folder
contains `metrics.json`, `summary.csv` and saved prediction profiles.
The metrics are global Pearson correlation, log-count Pearson correlation,
and mean Jensen–Shannon distance.

If a run fails, read its log and rerun the same command to reuse saved profiles.
For fresh predictions, change `full_run` to a new folder name. Do not reuse the
`quick_test` folder for a full run.

## Data and BED files

Signal and model files are downloaded into `data/atac/`:

- [ENCFF650PHA](https://www.encodeproject.org/files/ENCFF650PHA/): observed ATAC signal.
- [ENCFF113GSV](https://www.encodeproject.org/files/ENCFF113GSV/): ChromBPNet model weights.
- [ENCFF802KJC](https://www.encodeproject.org/files/ENCFF802KJC/): original IMR-90 peak regions.

The runner downloads the [shared hg38 FASTA](../../data/README.md#shared-human-reference) if needed.
