"""Prepare the reference-matched CAGI5 input and run configuration."""
from hashlib import md5, sha256
from pathlib import Path
import pandas as pd
import yaml

data = Path('data/mpra').resolve()
original = data / 'Kircher_GRCh38_ALL.tsv'
if md5(original.read_bytes()).hexdigest() != 'eacbb118660a3bd39575d88b1a03f7a7':
    raise ValueError('Unexpected CAGI5 input file')
source = pd.read_csv(original, sep='\t')
source.insert(0, 'source_row_id', source.index)
filtered = data / 'Kircher_GRCh38_REF_matched.tsv'
source.loc[source.source_row_id != 27570].to_csv(filtered, sep='\t', index=False)
config = yaml.safe_load(Path('benchmarks/mpra/cagi5.yaml').read_text())
config['data_dir'] = str(data)
out = Path('results/mpra').resolve()
out.mkdir(parents=True, exist_ok=True)
config['output_dir'] = str(out)
config['references']['hg38']['fasta_path'] = str(Path('data/genomes/hg38.fa').resolve())
settings = config['benchmarks']['cagi5_enformer_mapping']['config']
settings['variants_url'] = filtered.as_uri()
settings['variants_checksum'] = 'sha256:' + sha256(filtered.read_bytes()).hexdigest()
(out / 'cagi5.filtered.yaml').write_text(yaml.safe_dump(config))
