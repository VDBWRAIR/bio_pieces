import argparse
import subprocess
from glob import glob
import multiprocessing

from bio_bits.phylip import (
    make_renamed_phylip,
    rename_sequences,
)

def main():
    args = parse_args()
    # Script input options
    inputfasta = args.inputfasta
    phyml_options = args.phymloptions
    # Get renamed phylip and store the mapping of orig -> Seq#_ names
    mapping, phylipfile = make_renamed_phylip(inputfasta)
    # phyml options to run
    if '--no_memory_check' not in phyml_options:
        phyml_options += ' --no_memory_check -i {0}'.format(phylipfile)
    # Run phyml on renamed phylip
    run_phyml(phylipfile,phyml_options) 
    # Rename all phyml output files
    rename_list = glob(phylipfile + '_*')
    for f in rename_list:
        rename_sequences(f, mapping)

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'inputfasta',
        help='Fasta file to run through phyml(-i option)'
    )

    parser.add_argument(
        'phymloptions',
        default='',
        help='What options to pass on to phyml command.'
    )

    return parser.parse_args()

def run_phyml(phyml, phymloptions):
    '''
    Run phyml on a phymlfile
    '''
    # Run phyml
    phyml_cmd = 'phyml {0}'.format(phymloptions)
    print phyml_cmd
    try:
        subprocess.check_call(phyml_cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print '!!! There was an error running phyml !!!'
        print
        print e
        print
        print '----------------'

