import importlib.util
from pathlib import Path
import tempfile
import unittest
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
scorer = load('scorer', ROOT / 'benchmarks/expression/scripts/score_ct_specificity.py')

class ExpressionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        selected = self.base / 'selected.csv'
        pd.DataFrame({'gene_id':['ENSG00001.1,ENSG00002.1'], 'assembly':['human'],'on_targets':['ENCFFa'],'off_targets':['ENCFFb']}).to_csv(selected,index=False)
        self.selected = selected

    def test_scorer_aligns_shuffled_ids(self):
        truth = pd.DataFrame({'gene_id':['ENSG00001.1','ENSG00002.1'],'ENCFFa':[1.,4.],'ENCFFb':[3.,1.]})
        pred=truth.iloc[::-1].copy()
        score=scorer.score_predictions(truth,pred,self.selected,need_log=False)
        self.assertAlmostEqual(score['deviation_r'],1.)
    def test_scorer_uses_common_finite_rows(self):
        path = Path(self.selected)
        df = pd.read_csv(path);df.loc[0,'gene_id'] += ',ENSG00003.1';df.to_csv(path,index=False)
        truth = pd.DataFrame({'gene_id':['ENSG00001.1','ENSG00002.1','ENSG00003.1'],'ENCFFa':[1.,4.,np.nan],'ENCFFb':[3.,1.,2.]})
        pred = truth.fillna(3).iloc[::-1]
        score=scorer.score_predictions(truth,pred,path,need_log=False)
        self.assertEqual(score['genes_number'],2)
        self.assertAlmostEqual(score['deviation_r'],1.)

if __name__=='__main__': unittest.main()
