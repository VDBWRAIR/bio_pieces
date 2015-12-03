# -*- coding: utf-8 -*-

from future.utils import iteritems
from argparse import ArgumentParser
import os
import re
from Bio import SeqIO
from collections import defaultdict

def make_genemap(seqfile, delimiter, colnum, file_type):
    genemap = defaultdict(list)
    for record in SeqIO.parse( seqfile, file_type):
        parts = re.split(delimiter,record.description)
        gene = parts[colnum]
        genemap[gene].append( record ) 
    return genemap

def write_split_files( genemap, outformat, outdir):
    if not os.path.isdir(outdir):
        raise OSError("{0} is not a valid directory".format(outdir))
    for gene, records in iteritems(genemap):
        filename = "{0}/{1}.{2}".format(outdir, gene, outformat) 
        with open(filename, 'w') as outfile:
            SeqIO.write(records, outfile, outformat)

def parse_args():
    
    parser = ArgumentParser(
        description='Splits a fasta file into multiple files by ' \
            'grouping sequences by a specified column in the sequence ' \
            'identifier'
    )

    parser.add_argument(
        'seqfile',
        help='Sequence file path to split'
    )

    parser.add_argument(
        '--delimiter',
        '-d',
        default='__',
        help='The delimiter in the sequence id line to use'\
            '[Default: %(default)s]'
    )

    parser.add_argument(
        '--colnum',
        '-c',
        default=2,
        type=int,
        help='The column in the sequence identifier line to use to group, starting at 0'\
            '[Default: %(default)s]'
    )

    parser.add_argument(
        '--file-type',
        '-ift',
        default='fasta',
        choices=['fasta', 'fastq'],
        help='The file type the input file is[Default: %(default)s]'
    )
    parser.add_argument(
            '--out-format',
            choices=['fasta', 'fastq'],
            default=None
            )
    parser.add_argument( 
            '--outdir', '-o',
            required=True) 

    return parser.parse_args()

def split_and_write_files(seqfile, delimiter, colnum, outdir, file_type, out_format): 
    genemap = make_genemap( seqfile, delimiter, colnum, file_type)
    write_split_files( genemap, out_format, outdir) 

def main():
    args = parse_args()
    if args.out_format is None:
        args.out_format = args.file_type
    split_and_write_files(**args.__dict__ )


