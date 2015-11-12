#!/usr/bin/env python

##########################################################################
##                       Sequence Concat
##	Author: Tyghe Vallard						                        
##	Release Date: 5/30/2012							                            
##	Version: 1.1							                            
##	Description:							                            
##      This script will merge all sequences with identical names in a
##      file.
##
##  Example:
##      SampleFile1.fasta contents
##      >identifier1
##      AAAAAAAAAAAAAAAAAAAAAAAAA
##      >identifier1
##      TTTTTTTTTTTTTTTTTTTTTTTTT
##      >identifier1
##      GGGGGGGGGGGGGGGGGGGGGGGGG
##      >identifier2
##      CCCCCCCCCCCCCCCCCCCCCCCCC
##      >identifier2
##      AAAAAAAAAAAAAAAAAAAAAAAAA
##
##      OUTPUT
##      >identifier1
##      AAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGGGG
##      >identifier2
##      CCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAA
##      
##  VERSION HISTORY
##  -----------------
##  v1.1 - 6/08/2012
##    - Added options for doing gisaid formatted files
##    - Added support for different formatted required_segments
##########################################################################
import os
from optparse import OptionParser,OptionGroup
import sys
import cStringIO

from bio_bits import fasta
from bio_bits.fasta import UnknownIdentifierLineException

class SequenceConcat:
    _fasta_file_path = None
    _strip_chars = ""
    _parsed_fasta = None

    def __init__( self, fasta_file, file_type = 'genbank', strip_chars = "-" ):
        self._fasta_file_path = fasta_file
        self._strip_chars = strip_chars
        if file_type == 'genbank':
            try:
                self._read_genbank()
            except UnknownIdentifierLineException, e:
                print "An unknown identifier line was encountered in the fasta file. Is this a genbank file? If so use --type genbank"
                print e
                sys.exit( 1 )
        elif file_type == 'gisaid':
            try:
                self._read_gisaid()
            except UnknownIdentifierLineException, e:
                print "An unknown identifier line was encountered in the fasta file. Is this a gisaid file? If so use --type gisaid"
                print e
                sys.exit( 1 )

    def _read_genbank( self ):
        """ Reads the genbank file into a easy to work with dictionary """
        self._parsed_fasta = fasta.read_genbank_fasta( self._fasta_file_path, self._strip_chars )

    def _read_gisaid( self ):
        """ Reads the gisaid file into a easy to work with dictionary """
        self._parsed_fasta = fasta.read_gisaid_fasta( self._fasta_file_path, self._strip_chars )

    def prune_missing( self, fasta_dict, segments_expected = None ):
        """
            Prunes the dictionary of sequences so that only the sequences that have segments_expected
            amount of segments are included.

            Parameters:
                fasta_dict - Dictionary form of a fasta file from pyWRAIRLib.parser.fasta functions
                segments_expected - List of segments expected in what order(I.E [1,2,3,4,5,6,7,8] or ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'])
            
            Note:
                This is an inplace operation
            
            +++ Unit Test +++
            # Try an inhouse easy test
            >>> path = os.path.dirname( __file__ )
            >>> s = SequenceConcat( os.path.join( path, '../example_files/example1.txt' ), 'genbank' )
            >>> fasta = s.get_sequences()
            >>> pruned_fasta = s.prune_missing( fasta, range( 1, 4 ) )
            >>> print pruned_fasta
            {'ident4': {'1': 'AAAAAAAAAAAAAAAAA', '3': 'TTTTTTTTTTTTTTTTT'}, 'ident5': {'1': 'AAAAAAAAAAAAAAAAA'}, 'ident3': {'1': 'AAAAAAAAAAAAAAAAA', '2': 'CCCCCCCCCCCCCCCCC'}}
            >>> fasta != pruned_fasta
            True
            >>> len( fasta ) == 2
            True
            >>> s = SequenceConcat( os.path.join( path, '../example_files/example2.txt' ), 'gisaid' )
            >>> fasta = s.get_sequences()
            >>> pruned_fasta = s.prune_missing( fasta, ['A','B','C'] )
            >>> print pruned_fasta
            {'ident4': {'A': 'AAAAAAAAAAAAAAAAA', 'B': 'TTTTTTTTTTTTTTTTT'}, 'ident5': {'A': 'AAAAAAAAAAAAAAAAA'}, 'ident3': {'C': 'CCCCCCCCCCCCCCCCC', 'B': 'AAAAAAAAAAAAAAAAA'}}
            >>> fasta != pruned_fasta
            True
            >>> len( fasta ) == 2
            True
        """
        # We will return a dictionary that contains the sequences that have missing segments
        segments_missing = {}
        # Delete any of the sequences names from the dictionary that
        #  do not have the required amount of segments
        for seq_name in fasta_dict.keys():
            if len( fasta_dict[seq_name] ) != len( segments_expected ):
                # Copy the sequence
                segments_missing[seq_name] = fasta_dict[seq_name]
                # Delete it
                del fasta_dict[seq_name]
        return segments_missing

    def get_sequences( self ):
        """
            Return unmodified fasta dictionary
            
            +++ Unit Tests +++
            >>> s = SequenceConcat( '../example_files/example1.txt', 'genbank' )
            >>> len( s.get_sequences() ) > 1
            True
        """
        return self._parsed_fasta

    def get_merged_sequences( self, prune = True, segments_required = [1,2,3,4,5,6,7,8] ):
        """
            Returns a merged fasta formatted file
            Set prune to false if you don't want to prune out the sequences that don't have
            the correct amount of segments
        
            +++ Unit Tests +++
            >>> path = os.path.dirname( __file__ )
            >>> s = SequenceConcat( os.path.join( path, '../example_files/example1.txt' ), 'genbank' )
            >>> print s.get_merged_sequences( True, range( 1, 4 ) )
            >ident1
            AAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTT
            >ident2
            AAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC
            <BLANKLINE>
            >>> s = SequenceConcat( '../example_files/WR2848N.fasta', 'genbank' )
            >>> fh = open( '../example_files/WR2848N_merged.fasta' )
            >>> WR2848N_manually_merged = fh.read()
            >>> fh.close()
            >>> s.get_merged_sequences( ) == WR2848N_manually_merged[:-1]
            True
            >>> path = os.path.dirname( __file__ )
            >>> s = SequenceConcat( os.path.join( path, '../example_files/example2.txt' ), 'gisaid' )
            >>> print s.get_merged_sequences( True, ['A','B','C'] )
            >ident1
            AAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTT
            >ident2
            AAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC
            <BLANKLINE>
        """

        # The fasta sequences
        fasta_dict = self.get_sequences()

        # Fast String writer
        cs = cStringIO.StringIO()

        # Return variable
        output = None

        # Will hold missing sequences
        missing = None

        # If the prune option is set then prune the sequences
        if prune:
            # Make sure that the segments_required is in the right format
            # If string then split on ,
            if type( segments_required ) == str:
                segments_required = segments.required.split( ',' )
            # If already a list then it is ok
            elif type( segments_required ) == list:
                pass
            # Bail out if it gets here
            else:
                print( "Invalid value for required segments" )
                sys.exit( -1 )

            # Prune the dictionary
            missing = self.prune_missing( fasta_dict, segments_required )

            # Write a header to stderr
            if len( missing ):
                sys.stderr.write( "==================== Sequences Missing Segments  ====================\n" )

            #segments_required
            # Loop through the missing segment sequences
            for name, segments in missing.items():
                missing_segs = [str(i) for i in segments_required if str(i) not in segments]
                sys.stderr.write( ">%s is missing segment[s] %s\n" % (name, ",".join( missing_segs ) ) )

        # Loop through each sequence and merge the segments
        for name, segments in fasta_dict.items():
            cs.write( fasta.merge_segments( name, segments, segments_required ) )

        output = cs.getvalue()
        cs.close()
        
        return output
                
################ Script Functions #######################
def set_opts( parser ):
    """ Set script options """
    parser.add_option( "-f", "--fasta", dest="fasta_file", help="The fasta file of the sequences to merge" )
    parser.add_option( "-t", "--type", dest="db_type", default="genbank", help="What database type is this? gisaid and genbank are the only two options right now. Default: genbank" )
    parser.add_option( "--strip", dest="strip_chars", default="-", help="List of characters to strip from the sequences. Default is none" )
    parser.add_option( "--test", dest="test", action="store_true", help="Run the tests for this script" )
    prune_group = OptionGroup( parser, "Pruning Options" )
    prune_group.add_option( "--noprune", dest="prune", action="store_false", default=True, help="Don't prune out sequences that don't have the required amount of segments" )
    prune_group.add_option( "--segments_required", dest="segments_required", default="1,2,3,4,5,6,7,8", help="Required segments per sequence. See README for examples. Default: 1,2,3,4,5,6,7,8" )
    parser.add_option_group( prune_group )
    options,args = parser.parse_args()
    if not options.fasta_file and not options.test:
        parser.print_help()
        parser.error( "Need to specify the fasta file" )
    return options

# Run The script if this script is executed
def main():
    parser = OptionParser()
    options = set_opts( parser )
    if options.test:
        import doctest
        doctest.testmod()
    else:
        sc = SequenceConcat( options.fasta_file, options.db_type, options.strip_chars )
        print sc.get_merged_sequences( options.prune, options.segments_required.split( ',' ) )
