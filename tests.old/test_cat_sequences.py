from common import *

class Base(BaseTester):
    modulepath = 'bio_bits.cat_sequences'

    def setUp(self):
        super(Base, self).setUp()
        self.idlist = [
            '1__1_0',
            '1__2_1',
            '2__1_2',
            '2__2_3',
        ]
        self.seqrecs = []

    def _create_test_seqrecords(self, idlist):
        import random
        # Shuffle is in-place so copy first
        shuffled = [i for i in idlist]
        random.shuffle(shuffled)
        # Generate 10 random seqrecs
        seqrecs = []
        for id in shuffled:
            rec = self.make_seqrec(
                self.rand_seq(),
                self.rand_quals(),
                id
            )
            seqrecs.append(rec)
        return seqrecs

    def _check_sorted_seqs(self, seqrecs, idlist):
        for resultseq, expectid in zip(seqrecs, idlist):
            eq_( expectid, resultseq.id )

    def _create_mock_seqfile(self, seqrecs, filepath):
        from Bio import SeqIO
        with open(filepath,'w') as fh:
            SeqIO.write(seqrecs, fh, 'fasta')

class TestSplitSeqId(Base):
    functionname = 'split_seq_id'

    def test_splits_using_delimiter(self):
        id = [str(i) for i in range(5)]
        randseqrec = self.make_seqrec(
            self.rand_seq(),
            self.rand_quals(),
            '__'.join(id),
        )

        r = self._C(randseqrec, '__')
        eq_( id, r, "Did not recombine id" )

@attr('current')
class TestSortSequences(Base):
    functionname = 'sort_sequences'

    def setUp(self):
        super(self.__class__, self).setUp()
        self.seqrecs = self._create_test_seqrecords(self.idlist)

    def test_sorts_seqrecords(self):
        print self.seqrecs
        self._C(self.seqrecs, '__', [1,2])
        print self.seqrecs
        self._check_sorted_seqs(self.seqrecs, self.idlist)

    def test_sorts_using_single_item(self):
        print self.seqrecs
        self._C(self.seqrecs, '_(?<!_)', [1])
        print self.seqrecs
        self._check_sorted_seqs(self.seqrecs, self.idlist)

    def test_sortkey_outofbounds_raises_invalidsortkey(self):
        from bio_bits.cat_sequences import InvalidSortKey
        idlist = [
            '1__1__1',
            '1__1',
        ]
        seqrecs = self._create_test_seqrecords(idlist)
        try:
            self._C(seqrecs, '__', [3])
            ok_(False, 'Did not raise InvalidSortKey')
        except InvalidSortKey as e:
            eq_(
                '3 is an invalid sort key for the identifier 1__1',
                e.message
            )

    def test_sorts_case_insensitive(self):
        idlist = [
            'Apple_1',
            'Apple_2',
            'APPLE_3'
        ]
        seqrecs = self._create_test_seqrecords(idlist)
        self._C(seqrecs, '_', [1,2])
        self._check_sorted_seqs(seqrecs, idlist)

class TestSortSeqFiles(Base):
    functionname = 'sort_seq_files'

    def test_uses_basename_for_returned_dictionary_key(self):
        seqfiles = [abspath('seqfile1.fasta')]
        seqrecs = self._create_test_seqrecords(self.idlist)
        self._create_mock_seqfile(seqrecs, seqfiles[0])
        r = self._C(seqfiles, '__', [1,2])
        ok_( 'seqfile1.fasta' in r )
        self._check_sorted_seqs(r['seqfile1.fasta'], self.idlist)

class TestCatSeqRecords(Base):
    functionname = 'cat_seqrecords'

    def setUp(self):
        super(self.__class__,self).setUp()
        self.seqrecs = []
        self.seqrecs.append(self.make_seqrec(
            self.rand_seq(),
            self.rand_quals(),
            '1__2'
        ))
        self.seqrecs.append(self.make_seqrec(
            self.rand_seq(),
            self.rand_quals(),
            '1__1'
        ))

    def test_concats_sequences(self):
        r = self._C(self.seqrecs, '__', [1,2], True)
        expectseq = ''
        for rec in self.seqrecs:
            expectseq += str(rec.seq)
        eq_( expectseq, str(r.seq) )

    def test_creates_correct_sequence_id(self):
        r = self._C(self.seqrecs, '__', [1,2], True)
        eq_( '1', r.id )
        r = self._C(self.seqrecs, '__', [2,1], True)
        eq_( '1', r.id )

    def test_creates_correct_sequence_description_id(self):
        r = self._C(self.seqrecs, '__', [1,2], True)
        eq_( '1__2,1__1', r.description )

    def test_creates_correct_sequence_description_id_notkept(self):
        r = self._C(self.seqrecs, '__', [1,2], False)
        eq_( '', r.description )

class TestCombineSeqsInorder(Base):
    functionname = 'combine_seqs_inorder'

    @raises(Exception)
    def test_raises_error_if_all_seqs_not_same_length(self):
        seqrecs = self._create_test_seqrecords(self.idlist)
        self._create_mock_seqfile(seqrecs, 'seqf1.fasta')
        del seqrecs[-1]
        self._create_mock_seqfile(seqrecs, 'seqf2.fasta')

        seqfiles = ['seqf1.fasta','seqf2.fasta']
        seqperfile = self._run_func(
            self.modulepath,
            'sort_seq_files',
            seqfiles,
            '__',
            [1,2]
        )
        self._C(seqperfile, '__', [1,2])

    def test_combines_correctly(self):
        seqfiles = [
            'seqfile1.fasta',
            'seqfile2.fasta'
        ]
        seqrecs = self._create_test_seqrecords(self.idlist)
        for s in seqfiles:
            self._create_mock_seqfile(seqrecs, s)

        # Sort sequences now to compare later
        self._run_func(
            self.modulepath,
            'sort_sequences',
            seqrecs,   
            '__',
            [1,2]
        )

        seqperfile = self._run_func(
            self.modulepath,
            'sort_seq_files',
            seqfiles,
            '__',
            [1,2]
        )

        r = self._C(seqperfile, '__', [1,2])

        # zip together expected, result
        for eseqrec, rseqrec in zip(seqrecs,r):
            expectedseq = ''
            # Join together expected sequences(essentially double, triple...)
            for i in range(len(seqfiles)):
                expectedseq += str(eseqrec.seq)
            eq_( expectedseq, str(rseqrec.seq) )
