from common import *

class Base(BaseTester):
    modulepath = 'bio_bits.phyml_seqrename'

    def make_fasta_record(self, id, seq):
        fasta_fmt = '>{0}\n{1}\n'
        return fasta_fmt.format(id, seq)

    def writelines(self, lines, handle):
        ''' Writes lines to file and ensures newline '''
        if isinstance(handle, str):
            fh = open(handle, 'w')
        else:
            fh = handle
        for line in lines:
            if not line.endswith('\n'):
                line += '\n'
            fh.write(line)
        if isinstance(handle, str):
            fh.close()

class TestGetSeqEnumeration(Base):
    functionname = 'get_seq_enumeration'

    def test_gets_correct_enumeration(self):
        ids = ['Seq_{0}'.format(i) for i in range(10)]
        # Should be the resulting mapping
        enum = enumerate(ids)
        # Make mock fasta
        recs = [self.make_fasta_record(id,'ATGC') for id in ids]
        self.writelines(recs, 'input.fasta')

        r = self._C( 'input.fasta', 0 )
        for item in zip(enum, r):
            eq_( item[0][0], item[1][0] )
            eq_( item[0][1], item[1][1] )

class TestGetRenameList(Base):
    functionname = 'get_rename_list'

    def test_gets_same_order_as_input(self):
        # Reverse list
        ids = ['Seq_{0}'.format(i) for i in range(10)][::-1]
        self.writelines(ids, 'file')
        r = self._C('file', 'Seq_\d')
        eq_( ids, r )

class TestRenameList(Base):
    functionname = 'rename_list'

    def setUp(self):
        super(TestRenameList,self).setUp()
        self.rename_these = ['>Seq_{0}'.format(i) for i in range(10)][::-1]
        self.to_these = ['>Seq{0}'.format(i) for i in range(10)][::-1]
        self.renamefile = 'rename_these.txt'
        self.fromfile = 'from_these.txt'
        self.writelines(self.rename_these, self.renamefile)
        self.writelines(self.to_these, self.fromfile)

    def test_gets_rename_list(self):
        r = self._C(self.fromfile, self.renamefile, 'Seq_\d')
        
        e = zip(self.rename_these, self.to_these)
        for eitem, ritem in zip(e,r):
            # Compare without the >
            eq_( eitem[0][1:], ritem[0] )
            eq_( eitem[1][1:], ritem[1] )

    @raises(Exception)
    def test_exception_when_empty_renamelist(self):
        r = self._C(self.fromfile, self.renamefile, 'Seq_\d')

        self._C(self.fromfile, self.renamefile, 'Seq\d_')

    @raises(Exception)
    def test_exception_when_input_and_output_different_size(self):
        r = self._C(self.fromfile, self.renamefile, 'Seq_\d')
        with open(self.fromfile, 'a+') as fh:
            fh.write('>Seq_11\n')

        r = self._C(self.fromfile, self.renamefile, 'Seq_\d')

class TestRenameContents(Base):
    functionname = 'rename_contents'

    def test_renames_contents_of_file(self):
        rename_these = ['>Seq_{0}'.format(i) for i in range(10)][::-1]
        to_these = ['>Seq{0}'.format(i) for i in range(10)][::-1]
        self.writelines(rename_these, 'rename_these.txt')
        self.writelines(to_these, 'from_these.txt')

        from StringIO import StringIO
        outstream = StringIO()
        self._C(
            'from_these.txt',
            'rename_these.txt',
            'Seq_\d',
            outstream
        )

        contents = outstream.getvalue()
        for oldname in rename_these:
            ok_( oldname not in contents,
                'Found {0} still in file'.format(oldname)
            )

class TestMain(Base):
    functionname = 'main'

    def test_runs_main(self):
        inids = ['Seq_{0}'.format(i) for i in range(10)]
        recs = [self.make_fasta_record(id,'ATGC') for id in inids]
        self.writelines(recs, 'from.fasta')
        renameids = ['Seq{0}_'.format(i) for i in range(10)]
        recs = [self.make_fasta_record(id,'ATGC') for id in renameids]
        self.writelines(recs, 'rename.fasta')

        with patch('bio_bits.phyml_seqrename.argparse') as argparse:
            args = Mock(
                renamestr = 'Seq\d+_',
                fasta = 'from.fasta',
                renamefile = 'rename.fasta'
            )
            argparse.ArgumentParser.parse_args.return_value = args
            self._C()
