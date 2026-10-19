import tempfile,unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from prepare_evidence import calculate
from curation_export import export
from taxonomy_paths import resolve

class TestStages(unittest.TestCase):
    def test_robust_outlier(self):
        rows=[{'contig_id':f'c{i}','original_bin':'b1','gc':str(80 if i==9 else 40+i%3),'coverage':str(100 if i==9 else 10+i%3)} for i in range(10)]
        calculate(rows,'gc','coverage',None)
        self.assertEqual(rows[-1]['qc_outlier_flags'],'coverage,gc')
        self.assertFalse(rows[0]['qc_outlier_flags'])
    def test_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'asm.fasta').write_text('>A extra\nACGG\n>B\nGGTT\n')
            (p/'curated.tsv').write_text('contig_id\tcurated_bin\nA\tB1\nB\tB2\n')
            export(p/'asm.fasta',p/'curated.tsv',p/'out')
            self.assertIn('>A', (p/'out'/'bin_0001.fasta').read_text())
    def test_taxonomy(self):
        n,p,s=resolve('5',{'1':'1','5':'1'},{'1':'no rank','5':'species'},{'1':'root','5':'Example'})
        self.assertEqual(s,'resolved');self.assertIn('Example',p)

if __name__=='__main__':unittest.main()
