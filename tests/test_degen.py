from bio_bits.compat import unittest
import mock
from bio_bits.degen import Gene, get_gene_degen_overlap_info, main, csv_file_to_genes
from itertools import starmap
import sys
from nose.plugins.attrib import attr

genbank_id = 'KJ189367.1'
gb_file = 'tests/testinput/sequence.gb'
class DegenTest(unittest.TestCase):
    def setUp(self):
        self.seq = 'G'*83 + 'GGRAAY' + 420*'A' + 'RAGT'+'C'*8000 + 'YYY' +  'A' * 2000
        self.positions = [85, 88, 509, 8513, 8514, 8515]
        self.genestrs = ["anchored capsid protein"] * 2 + ["membrane glycoprotein precursor"] + ["nonstructural protein NS5"]*3
        self.nts = [self.seq[i] for i in self.positions]
        self.genbank_id = genbank_id
        self.expected_str =''
        for p, g, nt in zip(self.positions, self.genestrs, self.nts):
            self.expected_str += '\t'.join( map(str, [g, p, nt]) ) + '\n'
        self.expected_str = self.expected_str[:-1]
        self.geneinfo =    [ ["anchored capsid protein",             84, 426],
        ["membrane glycoprotein precursor", 426,923,     ],
        ["envelope protein",                924,2408,    ],
        ["nonstructural protein NS1",       2409,3464,   ],
        ["nonstructural protein NS2A",      3465,4118,   ],
        ["nonstructural protein NS2B",      4119,4508,   ],
        ["nonstructural protein NS3",       4509,6365,   ],
        ["nonstructural protein NS4A",      6366,6746,   ],
        ["2K peptide",                      6747,6815,   ],
        ["nonstructural protein NS4B",      6816,7562,  ],
        ["nonstructural protein NS5",        7563,10259,]]
        self.genes = list(starmap(Gene, self.geneinfo))


    #@mock.patch('Bio.SeqIO.parse')
    #NOTE: bio_bits/degen.py must exist or else Schema throws an error.
    @mock.patch('bio_bits.degen.parse_fasta')
    @mock.patch('sys.argv', ['_', 'bio_bits/degen.py', '--gb-id', genbank_id])
    @attr('download')
    def test_functional(self, mparse):
        '''Note: requires internet access to genbank.'''
        mparse.return_value = [mock.Mock(seq=self.seq)]
        #actual = next(get_result_table('_', ref_id=self.genbank_id))
        if not hasattr(sys.stdout, "getvalue"):
            self.fail("need to run in buffered mode")
        main()
        actual = sys.stdout.getvalue().strip() # because stdout is an StringIO instance
        self.assertMultiLineEqual(actual, self.expected_str)

    def test_get_gene_degen_overlap_info(self):
        expected = sorted(zip(self.genestrs, self.positions, self.nts))
        actual = sorted(get_gene_degen_overlap_info(self.genes, self.seq))
        self.assertListEqual(actual, expected)

    @mock.patch('bio_bits.degen.parse_fasta')
    #NOTE: bio_bits/degen.py must exist or else Schema throws an error.
    @mock.patch('sys.argv', ['_', 'bio_bits/degen.py', '--gb-file', gb_file])
    def test_functional_with_file(self, mparse):
        mparse.return_value = [mock.Mock(seq=self.seq)]
        if not hasattr(sys.stdout, "getvalue"):
            self.fail("need to run in buffered mode")
        main()
        actual = sys.stdout.getvalue().strip() # because stdout is an StringIO instance
        self.assertMultiLineEqual(actual, self.expected_str)

    def test_csv_to_gene(self):
        csvfile = 'tests/testinput/degen.csv'
        expected = [Gene(name='foo', start=1, end=2), Gene(name='bar', start=9, end=33)]
        actual = list(csv_file_to_genes(csvfile))
        self.assertListEqual(actual, expected)


    def test_tab_to_gene(self):
        csvfile = 'tests/testinput/degen.tab'
        expected = [Gene(name='foo', start=1, end=2), Gene(name='bar', start=9, end=33)]
        actual = list(csv_file_to_genes(csvfile))
        self.assertListEqual(actual, expected)



