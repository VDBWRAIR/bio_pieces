from StringIO import StringIO

from common import *
from mock import *
import pysam

from bio_bits import subsamplebam

class Base(BaseTester):
    modulepath = 'bio_bits.subsamplebam'

    def setUp(self):
        super(Base,self).setUp()

        self.refname = 'Ref1'
        self.reads = [
            'ATG',
            ' TGCA',
            'ATG',
            '  GCAT',
            '   CATG',
            '     TGC',
            '     TGC',
            'ATGCATGG'
        ]
        self.expected = [
            'Read1','Read2','Read6'
        ]
        self.mpysam = MagicMock()

        pileup_cols = self.build_mock_pileup(self.reads)
        # Mocked already called AlignmentFile
        self.mpysam.pileup.return_value = iter(pileup_cols)

        self.example_bam = join(EXAMPLES, 'example.bam')

    def build_mock_pileup(self, alignedreads):
        '''
        Build a mocked mpileup from a list of reads 

        # Here is an example of how you can build a mock pileup
        reads = [
            'ATG',
            ' TGCA',
            'ATG',
            '  GCAT',
            '   CATG',
            '     TGC',
            '     TGC',
            'ATGCATGC'
        ]
        x = build_mock_pileup(reads)
        
        Then x would be the same as if you had called pysam.AlignedFile.pileup()
        '''
        pileup_cols = []

        # Iterate across reference
        for i in range(len(alignedreads)):
            pile_col = Mock()
            pile_col.nsegments = 0
            pile_col.reference_id = ''
            pile_col.reference_pos = i
            pile_col.pileups = []

            # Build pileup for column
            # Iterate through each read
            for readnum, read in enumerate(alignedreads,start=1):
                # Only include reads where the read has a base at location
                if i >= len(read) or read[i] == ' ':
                    continue
                alignedsegment = Mock()
                alignedsegment.query_name = 'Read{0}'.format(readnum)
                alignedsegment.query_sequence = read.strip()

                pileup_read = Mock()
                pileup_read.alignment = alignedsegment
                pileup_read.query_position = i
                pile_col.pileups.append(pileup_read)
                pile_col.nsegments += 1
            pileup_cols.append(pile_col)

        return pileup_cols

    def mock_pysam_aligned_segment(self, *args, **kwargs):
        '''
        http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment

        Some common things to set
        reference_id
        reference_end
        reference_start
        reference_length
        query_sequence
        query_qualities
        query_name
        '''
        return Mock(**kwargs)

@attr('current')
class TestSubsetForReference(Base):
    functionname = 'get_readset_for_reference'

    def test_gets_correct_seqids(self):
        mpysam = MagicMock()
        pileup_cols = self.build_mock_pileup(self.reads)
        mpysam.pileup.return_value = iter(pileup_cols)

        with patch('bio_bits.subsamplebam.random') as mrandom:
            def random_sample(population, size):
                if size > len(population):
                    raise ValueError('size > population')
                return list(population)[:size]
            mrandom.sample = random_sample
            r = self._C(mpysam, 'reference', 2)
    
        expected = set(self.expected)
        eq_(expected, r)

    @raises(subsamplebam.MinimumDepthException)
    def test_position_not_enough_depth(self):
        mpysam = MagicMock()
        pileup_cols = self.build_mock_pileup(self.reads)
        mpysam.pileup.return_value = iter(pileup_cols)
    
        r = self._C(mpysam, 'reference', 10)

class TestFunctional(Base):
    @patch('bio_bits.subsamplebam.argparse')
    def test_writes_sam_to_output_stdout(self, margparse):
        args = Mock()
        args.input = self.example_bam
        args.samplesize = 10
        margparse.ArgumentParser.return_value.parse_args.return_value = args

        output = StringIO()
        with patch('bio_bits.subsamplebam.sys') as msys:
            msys.stdout = output
            subsamplebam.main()

        fh = open('out.sam','w')
        fh.write(output.getvalue())
        fh.close()
        # Make sure depth is correct across
        for pilecol in pysam.AlignmentFile('out.sam','r'):
            eq_(10, pilecol.nsegments)

    @patch('bio_bits.subsamplebam.argparse')
    def test_writes_sam_to_output_filepath(self, margparse):
        args = Mock()
        args.input = self.example_bam
        args.samplesize = 10
        args.output = 'test.out.sam'
        margparse.ArgumentParser.return_value.parse_args.return_value = args

        subsamplebam.main()

        fh = open('test.out.sam')
        samoutput = fh.read()
        fh.close()
        print samoutput
        # Make sure depth is correct across
        for pilecol in pysam.AlignmentFile('out.sam','r'):
            eq_(10, pilecol.nsegments)
