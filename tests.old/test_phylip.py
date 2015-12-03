from common import *

class Base(BaseTester):
    modulepath = 'bio_bits.phylip'

    def setUp(self):
        super(Base,self).setUp()

        # Number of sequences
        self.numseqs = 10

        # Orig -> Seq#_
        self.mapping = {}
        for i in range(self.numseqs):
            self.mapping[str(i)] = 'Seq{0}_'.format(i)

    def make_fasta(self, handle, seqids, sequences):
        ''' Write fasta sequences to handle '''
        for id,seq in zip(seqids,sequences):
            handle.write('>{0}\n{1}\n'.format(id,seq))

    def write_mapping_fasta(self, path, mapping, keyvalues='values'):
        # Same number of random seqs as generated sequences
        self.randseqs = [self.rand_seq() for i in range(self.numseqs)]

        # Sort the ids
        # This fetches either the keys(orig names) or values(new names)
        ids = sorted(getattr(mapping,keyvalues)())

        # Write the fasta
        with open(path,'w') as fh:
            self.make_fasta(fh, ids, self.randseqs)

        # Return the id, seq mapping
        return zip(ids, self.randseqs)


class TestRenameSequences(Base):
    functionname = 'rename_sequences'

    def test_replaces_contents_of_file(self):
        inputfasta = 'input.fasta'
        # Write fasta with all ids such as Seq#_
        idseq = self.write_mapping_fasta(inputfasta, self.mapping)
        # Get a mapping from that fasta
        idseqmap = dict(idseq)
        print idseqmap

        # Should rename all seq ids to 0...9
        r = self._C(inputfasta, self.mapping)

        with open(inputfasta) as fh:
            for line in fh: # Get id from id line
                if line.startswith('>'):
                    # Should be 0..9
                    id = line.rstrip()[1:]
                    # Gets the mapped Seq#_ name
                    mappedid = self.mapping[id]
                    # Get the sequence we would expect from before the rename
                    seq = idseqmap[mappedid]
                else: # Compare sequences
                    eq_( seq, line.rstrip(), '{0}: {1} != {2}'.format(id,seq,line.rstrip()) )

class TestMakeRenamedPhylip(Base):
    functionname = 'make_renamed_phylip'

    def test_mapping_is_correct(self):
        inputfasta = 'input.fasta'

        # Write a fasta with original names
        # Fasta is sorted by original names
        #   >0
        #   ATGC
        #   >1
        #   ATGC
        idseq = self.write_mapping_fasta(inputfasta, self.mapping, 'keys')
        # Run function
        # Mapping should be {'0':'Seq0_','1':'Seq1_',...}
        mapping, phyfile = self._C(inputfasta)

        # Mapping should be 0...9 and Seq0_...Seq9_
        for orig, new in mapping.items():
            # Original names are simply 0..9
            eq_( 'Seq{0}_'.format(orig), new )

    @raises(ValueError)
    def test_duplicate_seqids(self):
        inputfasta = 'input.fasta'
        with open(inputfasta,'w') as fh:
            fh.write('>1\nATGC\n')
            fh.write('>1\nATGC\n')
        self._C(inputfasta)

    def test_phylip_is_correct(self):
        inputfasta = 'input.fasta'

        # Write a fasta with original names
        idseq = self.write_mapping_fasta(inputfasta, self.mapping, 'keys')
        # Get a mapping from that fasta
        idseqmap = dict(idseq)

        # Run function
        mapping, phyfile = self._C(inputfasta)

        with open(phyfile) as fh:
            lineno = 0
            for line in fh:
                line = line.strip()
                if lineno == 0:
                    numseq, seqlen = line.split()
                    eq_( self.numseqs, int(numseq) )
                    eq_( 10, int(seqlen) )
                else:
                    id, seq = line.split()
                    # Original name is just the number portion of the name
                    orig = int(id.replace('Seq','').replace('_',''))
                    eq_( idseq[orig][1], seq )
                lineno += 1

class TestGetFastaSeqEnumerate(Base):
    functionname = 'get_fasta_seq_enumerate'

    def setUp(self):
        import random
        super(TestGetFastaSeqEnumerate,self).setUp()

        self.inputfasta = 'input.fasta'
        
        # Shuffle seq ids
        self.ids = [str(i) for i in range(10)]
        random.shuffle(self.ids)
        self.randseqs = [self.rand_seq() for i in range(self.numseqs)]
        # Write fasta with random order
        with open(self.inputfasta,'w') as fh:
            self.make_fasta(fh, self.ids, self.randseqs)

    def test_enumerates_from_start(self):
        count = 0
        for phyname, seqrec in self._C(self.inputfasta, start=0):
            expectedid = self.ids[count]
            expectedseq = self.randseqs[count]

            eq_( 'Seq{0}_'.format(count), phyname )
            eq_( expectedid, seqrec.id )
            eq_( expectedseq, str(seqrec.seq) )
            count += 1

    def test_enumerates_from_start_not_0(self):
        count = 10
        for phyname, seqrec in self._C(self.inputfasta, start=10):
            expectedid = self.ids[count-10]
            expectedseq = self.randseqs[count-10]

            eq_( 'Seq{0}_'.format(count), phyname )
            eq_( expectedid, seqrec.id )
            eq_( expectedseq, str(seqrec.seq) )
            count += 1

class TestGetSeqMapping(Base):
    functionname = 'get_seqmapping'

    def setUp(self):
        import random
        super(TestGetSeqMapping,self).setUp()

        self.inputfasta = 'input.fasta'
        
        # In order sequences 0...10
        self.ids = [str(i) for i in range(10)]
        self.randseqs = [self.rand_seq() for i in range(self.numseqs)]
        # Write fasta with random order
        with open(self.inputfasta,'w') as fh:
            self.make_fasta(fh, self.ids, self.randseqs)

        # Expected mapping for fastafile generated
        self.mapping = {}
        for id in self.ids:
            self.mapping['Seq{0}_'.format(id)] = id

    def test_gets_mapping_same_order_as_inputfile(self):
        r = self._C(self.inputfasta)

        eq_(len(self.mapping), len(r) )

        for expectnew, expectorig in self.mapping.items():
            eq_( r[expectnew], expectorig )

    @raises(ValueError)
    def test_duplicate_ids_raises_exception(self):
        inputfasta = 'input.fasta'
        with open(inputfasta,'w') as fh:
            fh.write('>1\nATGC\n')
            fh.write('>1\nATGC\n')
        self._C(inputfasta)
