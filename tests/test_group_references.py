import sys
import unittest
from bio_bits import group_references as pc
import io
from Bio import SeqIO
import os
from os.path import join, dirname, abspath, exists
import shutil
import mock
if sys.version[0] == '3':
    from io import StringIO as BytesIO
else:
    from io import BytesIO

THISD = dirname(abspath(__file__))

class TestParseContigs(unittest.TestCase):

    def setUp(self):
        self.columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
        self.samtext	=	'\n'.join([
            'read1\t1\tchr1\t1\t60	9M	=	1	1	TTTCGAATC	FFFFFFFFF	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq',
            'read2	1	*	1	60	8M	=	1	1	CTTCGATC	AFFDDDDD	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq',
            'read3	1	chr1	1	60	8M	=	1	1	CCGATCAA	FF@@@F@F	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq'
        ])
        self.reads = '''
@read1
TTTCGAATC
+
FFFFFFFFF
@read2
CTTCGATC
+
AFFDDDDD
@read3
CCGATCAA
+
FF@@@F@F
'''
        self.sam_extra_fields	=	'\n'.join([
            'read1\t1\tchr1\t1\t60	9M	=	1	1	TTTCGAATC	FFFFFFFFF	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq',
            'read2	1	chr2	1	60	8M	=	1	1	CTTCGATC	AFFDDDDD	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq\tfoo',
            'read3	1	chr1	1	60	8M	=	1	1	CCGATCAA	FF@@@F@F	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq\tfoo\tbar'
        ])
        self.seqrecs =  list(SeqIO.parse(BytesIO(self.reads), format='fastq'))
        self.sam_header = '\n'.join([
            '@HD	VN:1.3	SO:coordinate',
            '@SQ	SN:chr1	LN:10',
            '@SQ	SN:chr2	LN:10',
            '@RG	ID:MiSeq	SM:test	CN:None	PL:ILLUMINA',
        ])

    def test_sam_to_df(self):
        result = pc.samview_to_df(self.samtext)
        self.assertEquals(result.columns.tolist(), self.columns)
        self.assertEquals(result.ix[2]['QNAME'], 'read3')

    def test_sam_to_df_with_extra_fields(self):
        result = pc.samview_to_df(self.sam_extra_fields)
        self.assertEquals(result.columns.tolist(), self.columns)
        self.assertEquals(result.ix[2]['QNAME'], 'read3')
        self.assertEquals(result.ix[0]['QUAL'], 'FFFFFFFFF')
        self.assertEquals(result.ix[2]['QUAL'],   'FF@@@F@F')

    def test_main(self): #, margs):
        try:
            os.makedirs('tests/testoutput/pc')
        except:
            pass
        samfile = join(THISD, 'out.samtext')
        outdir = join(THISD, 'testoutput', 'pc')
        sys.argv = ['group_references', samfile, '--outdir', outdir]
        # Remove existing testdata if exists
        if exists(samfile):
            os.unlink(samfile)
        if exists(outdir):
            shutil.rmtree(outdir)
        with open(samfile, 'w') as out:
            out.write(self.sam_header + '\n')
            out.write(self.samtext)
        with mock.patch('bio_bits.group_references.sh') as mock_sh:
            mock_sh.samtools.return_value = self.samtext
            rcode = pc.main()
        self.assertEquals(0, rcode)
        expected_group1 = [self.seqrecs[0], self.seqrecs[2]]
        actual_group1 = SeqIO.parse(join(THISD, outdir, 'chr1.group.fq'), format='fastq')
        for le, ri in zip(expected_group1, actual_group1):
            self.assertEquals(*map(lambda r: (str(r.seq), r.letter_annotations, r.id), [le, ri]))
        actual_group2 = SeqIO.parse(join(THISD, outdir, 'unmapped.group.fq'), format='fastq')
        expected_group2 = [self.seqrecs[1]]
        for le, ri in zip(expected_group2, actual_group2):
            self.assertEquals(*map(lambda r: (str(r.seq), r.letter_annotations, r.id) , [le, ri]))
