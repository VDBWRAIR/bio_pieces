#!/usr/bin/env python2

import fileinput
import sys
import argparse

def get_csv_map(csvfile):
    '''
    Turn a 2-column csv into {key:value} dictionary

    '''
    mapping = {}
    with open(csvfile) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            l,r = line.rstrip().split(',')
            mapping[l] = r
    return mapping

def rename_fasta_identifiers(fastafiles, mapping, inplace=False):
    '''
    Take a list of fasta files or '-' and a mapping of {from1:to1, ..., fromN:toN}
    and rename all instances of keys on identifier lines to froms

    Setting inplace to True will edit the file in place otherwise output will go t
    standard out
    '''
    for line in fileinput.input(fastafiles, inplace=inplace):
        if not line.startswith('>'):
            # Just print sequence lines
            sys.stdout.write(line)
            continue
        # Split only one time
        p = line.rstrip().split(None, 1)
        # Id is first item before first space and without beginning >
        _id = p[0][1:]
        # In case the id is not in the mapping just replace it with itself
        if _id not in mapping:
            sys.stderr.write('{0} is not in provided mapping\n'.format(_id))
            newid = _id
        else:
            newid = mapping[_id]
        # Replace old id with new id
        sys.stdout.write(line.replace(_id, newid))

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'mappingfile',
        help='csv file with 2 columns From,To to use to rename fasta file ids'
    )
    
    parser.add_argument(
        'fastafile',
        default='-',
        help='Fasta file to rename or - to read from standard input'
    )

    parser.add_argument(
        '--inplace',
        default=False,
        action='store_true',
        help='Flag to determine if fastafile should be edited in place. Default ' \
            'is to output to standard output(terminal)'
    )

    return parser.parse_args()

def main():
    args = parse_args()
    mapping = get_csv_map(args.mappingfile)
    rename_fasta_identifiers(args.fastafile, mapping, args.inplace)
