from pathlib import Path
import unittest,json,hashlib,csv,ast,zipfile,xml.etree.ElementTree as ET
R=Path(__file__).resolve().parents[1]
class DeliveryTests(unittest.TestCase):
 def test_saved_results_lock(self):
  for n,h in json.loads((R/'results/SAVED_RESULTS_LOCK.json').read_text()).items():self.assertEqual(hashlib.sha256((R/'results'/n).read_bytes()).hexdigest(),h)
 def test_all_rows_and_intervals(self):
  rows=list(csv.DictReader((R/'results/COMPLETE_SENSITIVITY_TABLE.csv').read_text().splitlines()));self.assertEqual(len(rows),14)
  for r in rows:
   lo,hi=ast.literal_eval(r['CI95']);self.assertLess(lo,0);self.assertGreater(hi,0)
  primary=[r for r in rows if r['analysis']=='primary_30min'];self.assertEqual([(r['common_cases'],r['common_controls']) for r in primary],[('122','650'),('200','650')])
 def test_fixed_counts_and_precision(self):
  j=json.loads((R/'results/FIXED_PAIR_PRECISION_CHECK.json').read_text());self.assertFalse(j['material_problem']);self.assertEqual(j['unique_patients'],967);self.assertEqual(j['saved_outcome_reconstruction_max_error'],0);self.assertTrue(j['no_effect_function_or_bootstrap_called'])

if __name__=='__main__':unittest.main()
