import unittest
from os.path import dirname, basename, join, abspath
from nose.plugins.attrib import attr

import report

run = '/media/VD_Research/TMPDIR/RUNS/160204_M04171_0004_000000000-AJVMN'
samplename = '0126-2010'
_f = run + '/Data/Intensities/BaseCalls/{0}_S30_L001_R1_001.fastq.gz'.format(samplename)
_fr = 644208
_r = run + '/Data/Intensities/BaseCalls/{0}_S30_L001_R2_001.fastq.gz'.format(samplename)
_rr = _fr

class TestMiSeqRun(unittest.TestCase):
    def setUp(self):
        self.inst = report.MiSeqRun(run)
        
    def test_gets_fastq_gzs(self):
        sn = self.inst.get_samplenames()
        fqgz = self.inst.get_fastq_gz()
        self.assertEqual(2+2*len(sn), len(fqgz))

    def test_gets_paired_fastq_gzs(self):
        samplenames = self.inst.get_samplenames()
        fqgz = self.inst.get_fastq_gz()
        paired = self.inst.get_paired_fastq_gz()
        print paired
        self.assertEqual(len(fqgz)/2, len(paired))
        undet = paired[-1:]
        paired = paired[:-1]
        for forward, reverse in paired:
            self.assertIn(basename(forward).split('_')[0], samplenames)
            self.assertIn(basename(reverse).split('_')[0], samplenames)
        self.assertEqual(basename(undet[0][0]).split('_')[0], 'Undetermined')
        self.assertEqual(basename(undet[0][1]).split('_')[0], 'Undetermined')

    def test_gets_samplenames(self):
        self.assertTrue(self.inst.get_samplenames() > 2)

    def test_get_fqqs_per_samplename(self):
        x = self.inst.get_fqs_per_samplename()
        for sn, fr in x:
            f,r = fr
            self.assertTrue(f.basename().startswith(sn), "{0} does not start with {1}".format(f, sn))
            self.assertTrue(r.basename().startswith(sn), "{0} does not start with {1}".format(r, sn))

    @attr('slow')
    def test_read_stats(self):
        r,q,l,b = self.inst.read_stats(_f)
        assert _fr == r, "{0} != {1}".format(_fr, r)

    @attr('slow')
    def test_sample_read_stats(self):
        x = self.inst.sample_read_stats((_f,_r))
        r = x[0][0] + x[1][0]
        assert _fr+_rr >= r, "{0} != {1}".format(_fr+_rr, r)

    def test_sets_rundate(self):
        self.inst.rundate = '160204'

    @attr('slow')
    def test_integration_sample_read_stats(self):
        fqs = dict(self.inst.get_fqs_per_samplename())
        x = fqs[samplename]
        x = self.inst.sample_read_stats(x)
        r = x[0][0] + x[1][0]
        assert _fr+_rr >= r, "{0} != {1}".format(_fr+_rr, r)

    def test_run_stats(self):
        allstats = [
            ('s1', ((1,1,1,1),(1,1,1,1))),
            ('s2', ((2,2,2,2),(2,2,2,2))),
        ]
        x = self.inst.run_stats(allstats)
