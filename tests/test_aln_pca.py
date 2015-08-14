from __future__ import print_function

from StringIO import StringIO
from os.path import *

from . import unittest, THIS
import pandas as pd
import sh

from bio_pieces import aln_pca

aln = [
    ('>id1','ATGC'),
    ('>id2','AAGC'),
    ('>id3','CGTR'),
    ('>id4','ACGT')
]

id_matrix = pd.DataFrame(
    {
        'id1': [4.,3.,0.,2.],
        'id2': [3.,4.,0.,2.],
        'id3': [0.,0.,4.,0.],
        'id4': [2.,2.,0.,4.]
    },
    index=['id1','id2','id3','id4']
)

s_matrix = pd.DataFrame(
    {
        'A': [10,-8,-8,-8,0,1],
        'C': [-8,10,-8,-8,0,1],
        'G': [-8,-8,10,-8,0,1],
        'T': [-8,-8,-8,10,0,1],
        'R': [0,0,0,0,10,1],
        '~': [1,1,1,1,1,1],
    },
    index=['A','C','G','T','R','~']
)

class TestPairwiseIdentity(unittest.TestCase):
    def test_id_same(self):
        seq = 'A' * 10
        r = aln_pca.pairwise_identity(seq, seq)
        self.assertEqual(len(seq), r)

    def test_id_complete_diff(self):
        seq = 'A' * 10
        r = aln_pca.pairwise_identity(seq, seq.replace('A','C'))
        self.assertEqual(0, r)

    def test_id_half_diff(self):
        r = aln_pca.pairwise_identity('A'*10, 'A'*5+'C'*5)
        self.assertEqual(5, r)

    def test_diff_sizes(self):
        self.assertRaises(
            ValueError,
            aln_pca.pairwise_identity, 'A'*5, 'A'*4
        )

    def test_ignores_case(self):
        r = aln_pca.pairwise_identity('a'*10, 'A'*10)
        self.assertEqual(10, r)

    def test_utilizes_subst_matrix(self):
        r = aln_pca.pairwise_identity('ACGTR','AAAAA', s_matrix)
        e = s_matrix['A'].sum() - 1
        self.assertEqual(e, r)

    def test_subst_matrix_lookup_ignore_case(self):
        r = aln_pca.pairwise_identity('acgtr','AAAAA', s_matrix)
        e = s_matrix['A'].sum() - 1
        self.assertEqual(e, r)

    def test_substitution_matrix_missing_lookup_assigns_1(self):
        r = aln_pca.pairwise_identity('A'*5, 'A'*4+'-', s_matrix)
        self.assertEqual(10*4+1, r)

    def test_substitution_matrix_missing_allkey_raises_exception(self):
        smatrix = {
            'A': {'A': 1}
        }
        self.assertRaises(
            ValueError,
            aln_pca.pairwise_identity, 'A', '-', smatrix
        )

class TestIndexFasta(unittest.TestCase):
    def setUp(self):
        self.seq_string = '\n'.join(['\n'.join(x) for x in aln])
        self.aln_stringio = StringIO(self.seq_string)

    def test_ensure_ordered_index(self):
        r = aln_pca.index_fasta(self.aln_stringio)
        for _aln, _seq in zip(aln,r):
            search_id = _aln[0].replace('>','')
            self.assertEqual(_aln[1], r[search_id])

class TestIdentityMatrix(unittest.TestCase):
    def setUp(self):
        self.seq_string = '\n'.join(['\n'.join(x) for x in aln])
        self.aln_stringio = StringIO(self.seq_string)
        self.index_fasta = aln_pca.index_fasta(self.aln_stringio)
    
    def test_buils_correct_matrix(self):
        r = aln_pca.identity_matrix(self.index_fasta)
        self.assertTrue(id_matrix.equals(r))
        self.assertTrue(id_matrix.columns.identical(r.columns))
        self.assertTrue(id_matrix.index.identical(r.index))

class TestGeneratePlotFile(unittest.TestCase):
    def setUp(self):
        self.fasta = join(THIS, 'testinput', 'aln1.fasta')
        self.docdir = join(dirname(THIS), 'docs', '_static')

    def test_builds_simple(self):
        print(sh.aln_pca(self.fasta, outfile=join(self.docdir,'pca.png')))

    def test_accepts_subst_matrix(self):
        print(sh.aln_pca(
            self.fasta, s=join(dirname(THIS),'jalview_snm.txt'),
            outfile=join(self.docdir,'jalview.png')
        ))
