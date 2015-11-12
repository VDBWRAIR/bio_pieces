from os.path import *

from . import unittest, THIS
import mock

from bio_bits import beast_wrapper

class TestHoursStatesToSec(unittest.TestCase):
    def setUp(self):
        self.line = '20000   -61227.3011     -30325.1333     -30902.1678     932.046         1.60076E-4      0.35573         3487.00         3.39 hours/million states'
        self.chainlength = 500000000

    def test_converts_millions_to_sec(self):
        r = beast_wrapper.hours_states_to_sec(self.line, self.chainlength)
        self.assertEqual(6101755, r)

    def test_converts_billions_to_sec(self):
        line = self.line.replace('million','billion')
        r = beast_wrapper.hours_states_to_sec(line, self.chainlength)
        self.assertEqual(6101, r)

    def test_raises_value_error_for_invalid_beast_line(self):
        self.assertRaises(
            ValueError,
            beast_wrapper.hours_states_to_sec, 'foo', self.chainlength
        )

class TestSecToTime(unittest.TestCase):
    def test_converts_minutes(self):
        r = beast_wrapper.sec_to_time(1800)
        self.assertEqual('0d 00:30:00', r)

    def test_converts_hour(self):
        r = beast_wrapper.sec_to_time(3600)
        self.assertEqual('0d 01:00:00', r)

    def test_converts(self):
        r = beast_wrapper.sec_to_time(86400)
        self.assertEqual('1d 00:00:00', r)

    def test_converts_all(self):
        r = beast_wrapper.sec_to_time(86400+3600+1800+60+10)
        self.assertEqual('1d 01:31:10', r)

class TestGetXmlpathFromArgv(unittest.TestCase):
    def setUp(self):
        self.argv = ['beast', 'foo.xml', '-beagle_GPU']

    def test_gets_xmlpath(self):
        r = beast_wrapper.get_xmlpath_from_argv(self.argv)
        self.assertEqual('foo.xml', r)

    def test_gets_xmlpath_if_abspath(self):
        self.argv[1] = '/path/to/foo.xml'
        r = beast_wrapper.get_xmlpath_from_argv(self.argv)
        self.assertEqual('/path/to/foo.xml', r)

    def test_raises_value_error_if_missing_xmlpath(self):
        del self.argv[1]
        self.assertRaises(
            ValueError,
            beast_wrapper.get_xmlpath_from_argv, self.argv
        )

class TestGetChainglengthFromXml(unittest.TestCase):
    def setUp(self):
        self.xmlpath = join(THIS, 'beast.xml')
        self.xml = open(self.xmlpath)
        self.addCleanup(self.xml.close)

    def test_gets_chainlength(self):
        r = beast_wrapper.get_chainlength_from_xml(self.xml)
        self.assertEqual('1000000', r)

    def test_raises_value_error_if_missing_chainlength(self):
        self.xml = []
        self.assertRaises(
            ValueError,
            beast_wrapper.get_chainlength_from_xml, self.xml
        )

@mock.patch('bio_bits.beast_wrapper.sys.stdout.write')
class TestAddEstTimeToLine(unittest.TestCase):
    def setUp(self):
        self.lines = [
            'line1',
            '1 foo 1.0 hours/million states',
        ]

    def test_only_adds_sec_to_time_column_if_hours_per_line(self, mock_sout):
        beast_wrapper.add_est_time_to_line(1000000, self.lines[0])
        mock_sout.assert_called_once_with(
            'line1'
        )

    def test_adds_sec_to_time_column_if_hours_per_line(self, mock_sout):
        beast_wrapper.add_est_time_to_line(1000000, self.lines[1])
        mock_sout.assert_called_once_with(
            '1 foo 1.0 hours/million states\t0d 00:59:59\n'
        )

@mock.patch('bio_bits.beast_wrapper.sh')
class TestRunBeast(unittest.TestCase):
    def setUp(self):
        self.xmlpath = join(THIS, 'beast.xml')

    def test_calls_beast_correctly(self, mock_sh):
        argv = ['-beagle_SSE', self.xmlpath]
        func = mock.Mock()
        beast_wrapper.run_beast(*argv, _out=func)
        mock_sh.beast.assert_called_with(
            '-beagle_SSE',
            self.xmlpath,
            _out=func
        )

    def test_uses_sys_stdout_write_as_default_func(self, mock_sh):
        argv = ['-beagle_SSE', self.xmlpath]
        func = mock.Mock()
        with mock.patch('bio_bits.beast_wrapper.sys.stdout.write') as mock_sout:
            beast_wrapper.run_beast(*argv)
            mock_sh.beast.assert_called_with(
                '-beagle_SSE',
                self.xmlpath,
                _out=mock_sout
            )

@mock.patch('bio_bits.beast_wrapper.functools.partial')
@mock.patch('bio_bits.beast_wrapper.sys')
@mock.patch('bio_bits.beast_wrapper.sh')
class TestBeastWrapper(unittest.TestCase):
    def setUp(self):
        self.xmlpath = join(THIS, 'beast.xml')

    def test_outputs_modified_lines(self, mock_sh, mock_sys, mock_partial):
        mock_sys.argv = ('beast', '-beagle_SSE', self.xmlpath)
        beast_wrapper.beast_wrapper()
        args = mock_sh.beast.call_args
        self.assertEqual('-beagle_SSE', args[0][0])
        self.assertEqual(self.xmlpath, args[0][1])
        self.assertEqual(mock_partial.return_value, args[1]['_out'])
        mock_partial.assert_called_once_with(
            beast_wrapper.add_est_time_to_line, '1000000'
        )

class TestBeastEstTime(unittest.TestCase):
    def test_returns_correct_time(self):
        with mock.patch('bio_bits.beast_wrapper.sys') as mock_sys:
            mock_sys.argv.__getitem__.side_effect = [
                '1000000',
                ['1', 'foo', '1.0', 'hours/million', 'states']
            ]
            beast_wrapper.beast_est_time()
            mock_sys.stdout.write.assert_called_with('0d 00:59:59\n')
