'''

Usage:
    degen.py <fasta> (--gb-id <accession_id> | --gb-file <gbfile> | --tab-file <tabfile>)

Options:
    --gb-id=<accession_id>   Accession id for reference
    --gb-file=<gbfile>       Genbank file for reference
    --tab-file=<tabfile>     TSV/CSV file for reference with fields name,start,end
'''
from __future__ import print_function
from functools import partial
from collections import namedtuple
from itertools import starmap, product
from bio_pieces.compat import map, filter
import re
from bio_pieces.compat import StringIO
from Bio import Entrez, SeqIO

#for commandline stuff
from schema import Schema, Optional, Or
from docopt import docopt
import os
from operator import attrgetter as attr
import csv
from funcy.py3 import split
from funcy import compose

'''So GenBank can see how much you download.'''
Entrez.email = "micheal.panciera.work@gmail.com"


Gene = namedtuple('Gene', [ 'name', 'start', 'end'])

def seqrecord_to_genes(rec):
    '''
    :param Bio.SeqRecord rec: genbank record from SeqIO.parse format='genbank'
    :return iterable genes: iterator of gene objects (features with mat_peptied as their type)
    '''
    #Don't include `CDS`, that's whole-genome polypeptide
    GENE_TYPES = ('mat_peptide')
    genes = filter(lambda x: x.type in GENE_TYPES, rec.features)
    starts_ends_names = map(lambda f: ( f.qualifiers['product'][0], int(f.location.start), int(f.location.end), ), genes)
    return starmap(Gene, starts_ends_names)

def fetch_record_by_id(_id):
    return Entrez.efetch(db="nucleotide", id=_id, rettype='gb', retmode='text').read()

seq_parse_gb = partial(SeqIO.parse, format="genbank")
parse_fasta = partial(SeqIO.parse, format="fasta")
#assume genbank file only has one record (so use `next`)
id_to_record = compose(next, seq_parse_gb, StringIO, fetch_record_by_id)
id_to_genes = compose(seqrecord_to_genes, id_to_record)
genbank_file_to_genes = compose(seqrecord_to_genes, next, seq_parse_gb)
DEGENS = ['S', 'R', 'D', 'W', 'V', 'Y', 'H', 'K', 'B', 'M']
degen_positions = lambda seq:  (m.start() for m in re.finditer('|'.join(DEGENS), seq))

def row_to_gene(row):
    '''
    :param str row: row from csv.reader containing int1, int2, str (start, stop, name). header is ignored, and fields can be in any order,
    but the first integer must be less than the second.
    :return Gene gene: Gene object
    '''
    row = map(str.strip, row)
    digits, _gene_name  = split(str.isdigit, row)
    start, end = map(int, digits)
    assert start < end, "Start field should be first and less than end field. You supplied start %s end %s for gene %s" % (start, end, _gene_name[0])
    return Gene(next(_gene_name), start, end)

def open_generic_csv(csvfile):
    '''
    return an iterator (excluding header) of rows (as tuples) of a "csv" file with comma- or tab- seperators
    :param str csvfile: filename to csv/tsv file
    :return csv.Reader: iterator with header skipped
    '''
    dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters="\t,")
    csvfile.seek(0)
    reader = csv.reader(csvfile, dialect)
    has_header = csv.Sniffer().has_header(csvfile.read(1024))
    csvfile.seek(0)  # rewind
    if has_header: next(reader) #skip header
    return reader

csv_file_to_genes = compose(partial(map, row_to_gene), open_generic_csv, open)

def get_gene_degen_overlap_info(genes, seq):
    '''
    :param iterable genes: iterable of genes with attributes `start`, `end`, `name`
    :param str seq: nucleotide sequence
    :return list of tuples of form: (gene name, position, nt)... where degens and genes overlap.
    '''
    _degen_positions = degen_positions(str(seq))
    perms = product(genes, _degen_positions)
    result = ((gene.name, pos, seq[pos]) for gene, pos in perms if  gene.start <= pos <= gene.end)
    return result

def get_genes(ref_id=None, genbank_file=None, user_file=None):
    '''
    :param int ref_id: genbank accession id to get gene info from
    :param str genbank_file: filepath/filehandle for genbank file holding gene info
    :return iterable genes: iterable Gene objects with `start`, `end`, `name`
    '''
    assert  sum(map(bool, [ref_id, genbank_file, user_file])) == 1, "Must supply exactly one of accession id (%s) or gene_file (%s), or csv/tab-delimited file %s." % (ref_id, gene_file, tab_file)
    if ref_id:
        genes = id_to_genes(ref_id)
    elif genbank_file:
        genes = genbank_file_to_genes(genbank_file)
    elif user_file:
        genes = csv_file_to_genes(user_file)
    else:
        raise ValueError('Gene file or ref_id must be supplied.')
    return genes

def get_degen_gene_overlap(sequence, ref_accession_id=None, genbank_file=None, user_file=None):
    '''
    :param str sequence: dengue sequence
    :param int ref_accession_id: genbank accession id to get gene info from
    :param str genbank_file: filepath/filehandle for genbank file holding gene info
    :return iterable overlaps: iterable of tuples (`Gene name`, `Degen position`, `Degen base`)
    '''
    genes = get_genes(ref_accession_id, genbank_file)
    return get_gene_degen_overlap_info(genes, str(sequence))


'''
Functions for commandline app
'''
rowformat='{0}\t{1}\t{2}'
pretty_table = compose('\n'.join, partial(starmap, rowformat.format))

def main():
    scheme = Schema(
        { '<fasta>' : os.path.isfile,
         Optional('--gb-file') : Or(os.path.isfile, lambda x: x is None),
         Optional('--tab-file') : Or(os.path.isfile, lambda x: x is None),
         Optional('--gb-id') : Or(str, lambda x: x is None),
         })
    raw_args = docopt(__doc__, version='Version 1.0')
    args = scheme.validate(raw_args)
    fasta = parse_fasta(args['<fasta>'])
    genes = get_genes(args['--gb-id'], args['--gb-file'], args['--tab-file'])
    infos = map(partial(get_gene_degen_overlap_info, genes), map(attr('seq'), fasta))
    #need `list` to force evaluation of `print`
    list(map(print, map(pretty_table, infos)))

if __name__ == '__main__': main()
