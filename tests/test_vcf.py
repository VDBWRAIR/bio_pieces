import unittest
from bio_bits import vcfcat
import vcf
import sys
from bio_bits import vcfcat_main

class TestVcfCompare(unittest.TestCase):
    def setUp(self):
       self.infile = open('tests/testinput/947.vcf')
       self.smallfile = 'tests/testinput/small.vcf'
       self.ambiguousfile = 'tests/testinput/ambigous.vcf'
       self.recs = list(vcf.Reader(open(self.smallfile)))
       self.ambiguous_recs = list(vcf.Reader(open(self.ambiguousfile)))
       self.smalldifferent = 'tests/testinput/smalldifferent.vcf'
       self.diff_recs = list(vcf.Reader(open(self.smalldifferent)))

    def get_positions(self, result):
        return [rec.POS for rec in result]

#    def test_vcf_file_to_df(self):
#        result_df = vcfcat.vcf_file_to_df(vcf.Reader(self.infile))
#        raw_result_head = str(result_df.head()).strip()
#        result_head, expected_head = map(fix_string, [raw_result_head, self.expected_head])
#        self.assertEquals(expected_head, result_head)

    def test_filter_alt_equals(self):
        result = vcfcat.filter(self.recs, 'ALT', 'eq', ['C'])
        self.assertEquals(len(result) , 1)
        self.assertEquals(result[0].POS , 6)

    def test_filter_multi_alt_equals(self):
        result = vcfcat.filter(self.recs, 'ALT', 'eq', 'C, G')
        self.assertEquals(len(result) , 1)
        self.assertEquals(result[0].POS , 7)


    def test_filter_info_cbd(self):
        result = vcfcat.filter(self.recs, 'CBD', 'gt', 400)
        self.assertEquals(len(result) , 3)
        info = [(rec.POS, rec.REF) for rec in result]
        expected_info = [(8, 'A'), (9, 'G'), (10, 'T')]
        self.assertEquals(expected_info, info)


    def test_ambiguous_filter(self):
        result = vcfcat.ambiguous(self.ambiguous_recs)
        positions = self.get_positions(result)
        self.assertEquals([1, 8, 10], positions)
        cbs = [rec.INFO['CB'] for rec in result]
        self.assertEquals(['Y', 'R', 'Y'], cbs)


    def test_vcf_diff_str_no_threshold(self):
        result = vcfcat.diff(self.recs, self.diff_recs, 'REF', None)
        positions = self.get_positions(r for r in result)
        self.assertEquals([1, 10], positions)

    def test_vcf_diff_int_with_threshold(self):
        result = vcfcat.diff(self.recs, self.diff_recs, 'DP', 100)
        positions = self.get_positions(r for r in result)
        self.assertEquals([1, 3, 5], positions)

    def test_vcf_main(self):
        sys.argv = ['_', 'filter', self.infile.name, '--tag', 'CB', '--eq', 'A', '-c']
        result, _, _, __ = vcfcat_main.compute()
        self.assertEquals(3278, len(result))

    def test_exists(self):
        result = vcfcat.exists(self.recs, 'ALT')
        self.assertEquals(len(result), 2)
        existsin = 'tests/testinput/small.vcf'
        sys.argv = ['_', 'exists', existsin, '--tag', 'ALT', '--count']
        result, _, _, _,  = vcfcat_main.compute()
        self.assertEquals(len(result), 2)
        #TODO: could mock stdout here
        #vcfcat_main.main()
