from __future__ import print_function

import unittest
from bio_bits import sequence_split
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from  Bio.Alphabet import SingleLetterAlphabet
from collections import defaultdict
import mock
from argparse import Namespace

THISD = os.path.dirname(os.path.abspath(__file__))

class TestGeneMap(unittest.TestCase):
    def setUp(self):
        TESTFILES = os.path.join(THISD, 'expected')
        self.F_fastq = os.path.join(TESTFILES, 'F.fastq')
        self.example_fasta =  os.path.join(TESTFILES, 'example.fasta')
        self.gene1 =  os.path.join(TESTFILES, 'gene1.fasta')
        self.compile_genes = os.path.join(TESTFILES, 'compiled.fasta')
        self.F_fastq_expected_lengths = {"C1" : 23, "C2" : 25, "C3" : 27, "C4" : 50, "C5" : 75, "C6" : 26, "C7" : 12, "10" : 12}

    def test_genemap_fasta_gene1(self):
        record_1 = SeqRecord(seq=Seq('AAAA', SingleLetterAlphabet()), id='sample1__Brandomstuff__gene1__1', name='sample1__Brandomstuff__gene1__1', description='sample1__Brandomstuff__gene1__1', dbxrefs=[])
        record_2 = SeqRecord(seq=Seq('CCCC', SingleLetterAlphabet()), id='sample2__Arandomstuff__gene1__1', name='sample2__Arandomstuff__gene1__1', description='sample2__Arandomstuff__gene1__1', dbxrefs=[])
        expected =  {'gene1' : [record_1, record_2]}
        result = dict(sequence_split.make_genemap(self.gene1, '__', 2, 'fasta'))
        ''' BioSeq.Seq objects don't have proper __eq__ comparison implemented, so have to comapre __dict__ '''
        def seqs_equal(seq1, seq2):
            return str(seq1.seq) == str(seq2.seq) and seq1.id == seq2.id and seq1.name == seq2.name and seq1.description == seq2.description
        for key in expected.keys():
            self.assertTrue(key in result.keys())
            self.assertTrue(seqs_equal(expected[key][0], result[key][0]))

    def test_genemap_F_fastq_fields(self):
        result_map = sequence_split.make_genemap(self.F_fastq, ':', 1, 'fastq')
        expected_fields = ["@M02261:C1:000000000-A648D:1:1101:17702:1010 1:N:0:84", "NCTTGGCGTAAAGGCAGTGATTGCCGAGAGCTTTGAGCGAATACATCGTTCAAATTTAGTTGGTATGGGTATTTTACCGCTAACTTTNACCGGNNNNAATNCAAGATTAGATTTAAAGTTAGACGGTTCGGAAACTATTGATATTATAGGCTTAAGTGAACAAATAAAGCCTTATAATCCTGTTAAATGCATGATAA", "+", "#8ACCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG#:DFGG####::D#:CFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"]
        result_fields = result_map['C1'][0].format("fastq").strip().split('\n')
        self.assertEquals(expected_fields, result_fields)

    def test_genemap_F_fastq_compiled(self):
       result_map = sequence_split.make_genemap(self.F_fastq, ':', 1,'fastq')
       result_lengths = dict( (seqname, len(seqs)) for seqname, seqs in result_map.items())
       self.assertEquals(8, len(result_map))
       self.assertEquals(self.F_fastq_expected_lengths, result_lengths)

@mock.patch('bio_bits.sequence_split.ArgumentParser.parse_args')
class TestFunctional(unittest.TestCase):

    def setUp(self):
        self.outdir =  os.path.join(THISD, 'testoutput')
        try:
            os.mkdir(self.outdir)
        except:
            pass
        self.F_fastq_expected_lengths = {"C1" : 23, "C2" : 25, "C3" : 27, "C4" : 50, "C5" : 75, "C6" : 26, "C7" : 12, "10" : 12}
        self.default_args = Namespace(seqfile='tests/F.fastq', delimiter=':', file_type='fastq', outdir=self.outdir, colnum=1, out_format=None)

    def yield_files_and_lengths(self, outformat):
        for seq_name, num_seqs in self.F_fastq_expected_lengths.items():
            filename = os.path.join(self.outdir, "{0}.{1}".format(seq_name, outformat))
            yield filename, num_seqs

    def test_F_fastq_compiled_files_exist(self, m_parse_args):
        m_parse_args.return_value = self.default_args
        sequence_split.main()
        for filename, _ in self.yield_files_and_lengths('fastq'):
            self.assertTrue(os.path.isfile(filename))

    def test_F_fastq_to_fasta_files_exist(self, m_parse_args):
        m_parse_args.return_value =  Namespace(seqfile='tests/F.fastq', delimiter=':', file_type='fastq', outdir=self.outdir, out_format='fasta', colnum=1)
        sequence_split.main()
        for filename, _ in self.yield_files_and_lengths('fasta'):
            self.assertTrue(os.path.isfile(filename))

    def test_F_fastq_to_fasta_num_records(self, m_parse_args):
        m_parse_args.return_value =  Namespace(seqfile='tests/F.fastq', delimiter=':', file_type='fastq', outdir=self.outdir, out_format='fasta', colnum=1)
        sequence_split.main()
        for filename, num_records in self.yield_files_and_lengths('fasta'):
            result_records = list(SeqIO.parse(filename, 'fasta'))
            self.assertEquals(num_records, len(result_records))

    def test_F_fastq_to_fasta_some_lines_exist(self, m_parse_args):
        m_parse_args.return_value =  Namespace(seqfile='tests/F.fastq', delimiter=':', file_type='fastq', outdir=self.outdir, out_format='fasta', colnum=1)
        sequence_split.main()
        expected_lines = [">M02261:C1:000000000-A648D:1:1101:17702:1010 1:N:0:84",
"NCTTGGCGTAAAGGCAGTGATTGCCGAGAGCTTTGAGCGAATACATCGTTCAAATTTAGTTGGTATGGGTATTTTACCGCTAACTTTNACCGGNNNNAATNCAAGATTAGATTTAAAGTTAGACGGTTCGGAAACTATTGATATTATAGGCTTAAGTGAACAAATAAAGCCTTATAATCCTGTTAAATGCATGATAA", ">M02261:C1:000000000-A648D:1:1101:15313:1016 1:N:0:84"," NAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAANAGAAGNNAANAGAAGAGAAGAGAAGAGAGAAGAGGAGAGACA"]
        expected = ''.join(expected_lines).replace(' ', '')
        self.assertTrue(expected in open(os.path.join(self.outdir, 'C1.fasta')).read().replace('\n', '').replace(' ', ''))


    @mock.patch('bio_bits.sequence_split.os.path.isdir', return_value=False)
    def test_outdir_exists_raises_error(self, m_isdir, m_parse_args):
        m_parse_args.return_value = self.default_args
        self.assertRaises(OSError, sequence_split.main)

    def test_F_fastq_compiled_correct_record_number(self, m_parse_args):
        m_parse_args.return_value = self.default_args
        sequence_split.main()
        for filename, expected_num_records in self.yield_files_and_lengths('fastq'):
            result_records = list(SeqIO.parse(filename, 'fastq'))
            self.assertEquals(expected_num_records, len(result_records) )













