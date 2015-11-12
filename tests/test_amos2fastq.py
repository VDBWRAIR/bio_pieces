try:
    import unittest2 as unittest
except ImportError:
    import unittest

from operator import itemgetter
from pandas.util.testing import assert_frame_equal
import pandas as pd
import mock
from bio_bits import   amos2fastq_main
from bio_bits import amos2fastq as a2f
import os
import random
from functools import partial
''' Python3 compatibility '''
from past.builtins import map, filter

#TODO: Restore a second input file.
#TODO: Figure out a better way to produce test inputs. SeqIO?
'''
Assumptions:
    1. Each RED corresponds to exactly one FASTQ entry.
Want to test that we handle cases like:
    1. The same read maps to two contigs
    2. Two REDS/fastq entries have the same sequence string
'''
THISD = os.path.dirname(os.path.abspath(__file__))
inputdir = os.path.join(THISD, 'testinput')
outputdir = os.path.join(THISD, '.')
expecteddir  = os.path.join(THISD, 'expected')

def make_seq_record(seq):
    m = mock.Mock()
    m.seq.tostring.return_value = seq
    m.some_id = random.randint(0, 999999)
    return m

def make_mock_red(iid, eid, seq, qlt):
    return mock.Mock(**locals())

class TestAmos2Fastq(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_seqs_by_ctg(self):
        iids_by_ctg = [[1, 2, 4], [1, 2], [3, 5], [4], [6], [2, 4]]
        get_seq = lambda i: 'ACGT' if i in [2, 4] else 'A'*i
        reds = [make_mock_red(i, i, get_seq(i), 'D'*len(get_seq(i)))
                         for i in range(1, len(iids_by_ctg)+ 1)]
        fastq_records = [make_seq_record(get_seq(i)) for i in range(1, len(iids_by_ctg)+1)]
        saved_records = fastq_records[:]
        random.shuffle(fastq_records)
        result = filter(bool, map(tuple, a2f.get_seqs_by_ctg(fastq_records, reds, iids_by_ctg)))
        getters_by_ctg = map(lambda A: itemgetter(*[i - 1 for i in A]), iids_by_ctg)
        raw_expected = [get(saved_records) for get in getters_by_ctg]
        expected = filter(bool, map(lambda a: (a,) if type(a) is not tuple else a, raw_expected))
        # Bio.Seq objects can't be compared directly need to get __dict__ attribute
        dicter = lambda tup: [obj.seq.tostring() for obj in tup]
        self.assertEquals(map(dicter, expected), map(dicter, result))
        pass

class TestAmos2FastqMain(unittest.TestCase):

    def setUp(self):
        input_join = partial(os.path.join, inputdir)
        expected_join = partial(os.path.join, expecteddir)
        outnames = ['contig-0.fastq', 'contig-1000000.fastq', 'contig-2000000.fastq']
        self.fastqs = map(input_join, ['1.R1.unmap.fastq', '1.R2.unmap.fastq'])
        self.amos = input_join('regression.afg')
        self.expected_files = map(expected_join, outnames)
        self.result_files = outnames #map(output_join, outnames)

    @mock.patch('bio_bits.amos2fastq_main.docopt')
    def test_main(self, mockopt):
        mockopt.return_value = {'<fastqs>' : self.fastqs, '--amos' : self.amos}
        read_file = lambda fn, fn2: (open(fn).read(), open(fn2).read())
        amos2fastq_main.main()
        map(self.assertTrue, map(os.path.exists, self.result_files))


def make_df_index0(rows):
    return pd.DataFrame(rows)

class TestPandasUtils(unittest.TestCase):

    def test_df_from_collection_attributes(self):
        import mock
        mocks = [mock.Mock() for i in range(5)]
        [[mock_do() for i in range(index)] for index, mock_do in enumerate(mocks)]
        columns = ['call_count', 'called']
        expected = pd.DataFrame( [(i, bool(i)) for i in range(5)])
        expected.columns = columns
        result = a2f.df_from_collection_attributes(columns, mocks)
        assert_frame_equal(expected, result)

    def test_series_contains_nan(self):
        has_nan = pd.Series([0, 1, 2, float('nan'), 3, 4])
        no_nan = pd.Series([2, 4, 5, 'a'])
        self.assertTrue(a2f.series_contains_nan(has_nan))
        self.assertFalse(a2f.series_contains_nan(no_nan))

#    def test_join_non_unique_dataframes(self):
#        '''
#        df1 and df2 share an index with duplicates, check that it is aligned correctly
#        '''
#        rows1 = [('A', 'A1'), ('B', 'B1'), ('A', 'A2'), ('C', 'C1')]
#        rows2 = [('A', '0A', False), ('B', '0B', True), ('A', '00A', False), ('C', '00C', True)]
#        self.df1, self.df2 = map(make_df_index0, (rows1, rows2))
#        self.df1.columns = ['0', '1']
#        self.df2.columns = ['0', '1', '2']
#        self.df1, self.df2 = self.df1.set_index('0'), self.df2.set_index('0')
#        result = a2f.join_non_unique_dataframes(self.df1, self.df2)
#        expected = pd.DataFrame(
#                [('A', 0, 'A1', '0A', True), ('B', 0, 'B1', '0B', True),
#                ('A', 1, 'A2', '00A', False), ('C', 0, 'C1', '00C', True)]
#                ).set_index(0).set_index(1, append=True)
#        assert_frame_equal(result, expected)

