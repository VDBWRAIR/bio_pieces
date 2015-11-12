from __future__ import print_function
from . import unittest, THIS
from os.path import *
import sys
import mock
import sh

import pandas as pd

from bio_bits import beast_checkpoint

class Base(unittest.TestCase):
    def setUp(self):
        self.testlog = join(THIS,'beast.log')
        self.testxml = join(THIS,'beast.xml')
        self.testtrees = join(THIS,'beast.trees')

class BasePandas(Base):
    def setUp(self):
        super(BasePandas,self).setUp()
        self.patch_pandas = mock.patch('bio_bits.beast_checkpoint.pd')
        self.mock_pandas = self.patch_pandas.start()
        self.addCleanup(self.patch_pandas.stop)

class TestParseLogfile(BasePandas):
    def test_calls_pandas_correctly(self):
        r = beast_checkpoint.parse_logfile('foo.log')
        args = self.mock_pandas.read_csv.call_args
        self.assertEqual('foo.log', args[0][0])

    def test_returns_pandas_dataframe(self):
        self.mock_pandas.read_csv.return_value = 'bar'
        r = beast_checkpoint.parse_logfile('foo.log')
        self.assertEqual('bar', r)

class TestHashLastLogEntries(BasePandas):
    def test_sets_all_param_values_from_last_entry(self):
        df = pd.read_csv(self.testlog, sep='\t', comment='#')
        self.mock_pandas.read_csv.return_value = df
        r = beast_checkpoint.hash_last_log_entries('foo.log','bar.log')
        self.assertEqual(30000, r['state'])
        self.assertEqual(88.06244116021935, r['speciation'])

    def test_returns_same_order_as_log_header(self):
        df = pd.read_csv(self.testlog, sep='\t', comment='#')
        self.mock_pandas.read_csv.return_value = df
        r = beast_checkpoint.hash_last_log_entries('foo.log','bar.log')
        f = list(filter(lambda k: k.startswith('freq'), r.keys()))
        self.assertEqual(['frequencies1','frequencies2','frequencies3','frequencies4'], f)

class TestFunctional(Base):
    def test_mods_parameter_id_lines(self):
        log = pd.read_csv(self.testlog, sep='\t', comment='#')
        cmd = '{script} {xml} {trees} {log}'.format(
            script=join(dirname(THIS),'bio_bits','beast_checkpoint.py'),
            xml=self.testxml, trees=self.testtrees, log=self.testlog
        )
        ac = log.tail(1)['ac'].values[0]
        freqs = log.tail(1)[
            ['frequencies1','frequencies2','frequencies3','frequencies4']
        ]
        freqs = [str(v.values[0]) for k,v in freqs.iteritems()]
        freqs = ' '.join(freqs)

        # These all need to exist in output
        find_in_lines = [
            '<parameter id="frequencies" value="0.225554'.format(freqs),
            '<parameter id="ac" value="{0}" lower="0.0" upper="Infinity"/>'.format(ac),
        ]

        not_find_in_lines = [
            '<parameter id="treeModel.rootHeight" value=',
        ]

        found = []
        def find_lines(line):
            ''' look for required lines in output '''
            for l in find_in_lines + not_find_in_lines:
                if l in line:
                    found.append(l)

        out = sh.python(cmd.split(), _out=find_lines, _err=find_lines)
        print("Looking for lines {0}".format(find_in_lines))
        print("Found {0}".format(found))
        self.assertEqual(sorted(find_in_lines), sorted(found))
