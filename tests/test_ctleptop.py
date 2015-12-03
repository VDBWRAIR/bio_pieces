from bio_bits.compat import unittest, StringIO
from bio_bits import ctleptop
import mock
import os
from bio_bits.ctleptop import create_args
#from sys import version_info
# if version_info.major == 2:
# import __builtin__ as builtins
# else:
# import builtins

THIS = os.path.dirname(os.path.abspath(__file__))


class TestCtleptopFunctionTest(unittest.TestCase):
    infile = os.path.join(THIS, "Den4_MAAPS_TestData16.fasta")
    outfile = "out_file.fa"

    """Base TestCase class, set up a CLI parser"""

    def setUp(self):
        #self.patch_fileinput = mock.patch('bio_bits.ctleptop.access_mixed_aa')
        #self.mock_fileinput = self.patch_fileinput.start()
        #self.mock_fileinput_input = self.mock_fileinput.input
        # self.addCleanup(self.mock_fileinput.stop)
        #self.patch_sys_stdout = mock.patch('bio_bits.ctleptop.sys.stdout')
        #self.patch_sys_stderr = mock.patch('bio_bits.ctleptop.sys.stderr')
        #self.mock_stdout = self.patch_sys_stdout.start()
        #self.mock_stderr = self.patch_sys_stderr.start()
        #self.addCleanup(self.mock_stdout.stop)
        #self.addCleanup(self.mock_stderr.stop)
        #self.patcher_argparse = mock.patch('bio_bits.ctleptop.create_args')
        #self.mock_argparse = self.patcher_argparse.start()
        # self.addCleanup(self.mock_argparse.stop)
        #parser = ctleptop.create_args()
        #self.parser = parser
        pass

    #@mock.patch('bio_bits.ctleptop.sys.argv')
    #import sys
    @mock.patch('sys.argv', ['DUMMY', '-i', infile, '-o', outfile])
    def test_args(self):
        # ['-i' , self.infile, '-o', self.outfile])
        parser = ctleptop.create_args()
        #self.parser = parser
        self.assertEquals(parser.i, self.infile)
        self.assertEquals(parser.o, self.outfile)
'''
    @mock.patch('bio_bits.ctleptop.create_args')
    #@mock.patch('bio_bits.ctleptop.open_f')
    def test_main(self, call_main):
        from StringIO import StringIO
        args = ctleptop.create_args()
        #outfile = args.o
        ctleptop.main()
        try:
            with mock.patch('__builtin__.open', return_value = StringIO(args.o)):
                fileout =args.o
                self.assertEqual(fileout, args.o)
        except ImportError:
            print "Unable to open outfile"
        self.assertEquals(args.i, self.infile)
        #call_main.assert_called_once_with(args.o)
    @mock.patch('Bio.SeqIO.parse')
    def test_access_mixed_aa(self, mock_mixed_aa):
        for rec in mock_mixed_aa(self.infile, 'fasta'):
            header, seqline = rec.id, str(rec.seq)
            self.assertTrue(isinstance(">", header))
            self.assertTrue(isinstance("AA", seqline))
            mock_mixed_aa.assert_called_once_with(self.infile, 'fasta')
        ctleptop.access_mixed_aa(self.infile)
'''
class TestSetOptionDebug(unittest.TestCase):
    def setUp(self):
        pass
        @mock.patch('sys.stdout', new_callable=StringIO)
        def test_options_debug(self, stdout_mock):
            val = "Start processing and writing the output file to please please wait ... "
            ctleptop.main()
            self.assertEquals(stdout_mock.getvalue(), val)

class TestCtleptopTest(unittest.TestCase):

    """Docstring for CtleptopTest. """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getaalist(self):
        result = ctleptop.getaalist(['AAA', 'ACT'])
        self.assertEquals(result, ['K', 'T'])

    def test_nearbyPermutations(self):
        result = ctleptop.nearbyPermutations("AAR")
        self.assertEquals(result, set(['AAG', 'AAA']))
        result2 = ctleptop.nearbyPermutations("ARR")
        self.assertEquals(result2, set(['AGG', 'AAG', 'AAA', 'AGA']))
        result3 = ctleptop.nearbyPermutations("AAA")
        self.assertEquals(result3, set(['AAA']))

    def test_permutations(self):
        """Test listReplace function.
        :returns: bool

        """
        result = ctleptop.permutations(set(['CA']), ['A', 'T'])
        self.assertEquals(result, set(['ACA', 'TCA']))

    def test_list_overlap(self):
        """Check if two list overlap
        :returns: bool

        """
        result = ctleptop.list_overlap(
            "GAY", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(True,  result)
        result2 = ctleptop.list_overlap(
            "WTA", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(True,  result2)
        result3 = ctleptop.list_overlap(
            "TTT", ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M'])
        self.assertEquals(False,  result3)

    def test_isGap(self):
        result = ctleptop.isGap(["K/R", "I/T"], ["ARG", "NNN"])
        self.assertEquals(result, ["K/R", "GAPFOUND"])
