'''
Usage: group_references <samfile> [--outdir <DIR>]

Options:
    --outdir=<DIR>,-o=<DIR>   output directory [Default: group_references_out]

Create separate fastq files for each reference in a samfile.
'''

from docopt import docopt
import pandas as pd
from schema import Schema, Use
import sys
import os
import sh
if sys.version[0] == '3':
    from io import StringIO as BytesIO
else:
    from io import BytesIO
import string
import re
sam_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
def samview_to_df(rawtext):
    '''
    :param str rawtext: text of .sam file or output of `samtools view`
    :return: pandas.DataFrame 
    '''
    samtext = '\n'.join( fixline(row) for row in str(rawtext).split('\n') )
    as_bytes = BytesIO(samtext)
    #Write test to ensure handles varaibles #columns
    return pd.read_csv(as_bytes, names=sam_columns, usecols=sam_columns, delimiter='\t', header=None, squeeze=True)

def fixline(row):
    '''
    Fix lines of `samtools view` which are unmapped so they can be easily parsed by pandas.
    :row str row: a line from `samtools view`
    :return: fixed line
    '''
    newrow = []
    cols = row.split('\t')[:len(sam_columns)]
    if len(cols) > 1 and cols[2] == '*':
        cols[2] = 'unmapped'
    return '\t'.join(cols)

def get_seqs_by_ctg(outdir, rawtext):
    '''
    Splits a given sam view into seperate fastq files, grouped by references (the first column RNAME), and saves to outdir.
    :param str outdir: directory result fastqs are saved to
    :param str rawtext: output of `samtools view`
    :return: None
    ''' 
    sam_df = samview_to_df(rawtext)
    contig_groups = sam_df.groupby('RNAME')
    fastq = "@{0}\n{1}\n+\n{2}".format
    for group in contig_groups:
        _ref, reads = group[0], group[1]
        ref = re.sub('[%s]' % string.punctuation, '_', _ref) 
        with open("{0}/{1}.group.fq".format(outdir, ref), 'w') as out:
            map(out.writelines, '\n'.join(map(fastq, reads.QNAME, reads.SEQ, reads.QUAL)))
            out.write('\n')

def main():
    '''
    Call `samtools view` on the input file and split into fastqs by RNAME column.
    '''
    raw_args = docopt(__doc__)
    scheme = Schema({
        '<samfile>' : str,
        '--outdir' : str})
    parsed_args = scheme.validate(raw_args)
    outdir = parsed_args['--outdir']
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    infile = parsed_args['<samfile>']
    view = str(sh.samtools('view', infile, S=True)) if infile.endswith('.sam') else str(sh.samtools('view', infile))
    get_seqs_by_ctg(outdir, view)
    return 0

if __name__ == '__main__':
    main()
