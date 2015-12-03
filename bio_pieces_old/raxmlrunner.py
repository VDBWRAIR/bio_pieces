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
    outputsuffix = args.outputsuffix
    raxml_options = args.raxmloptions
    # Get renamed phylip and store the mapping of orig -> Seq#_ names
    mapping, phylipfile = make_renamed_phylip(inputfasta)
    # Raxml options to run
    raxml_options += ' -n {0} -s {1} --no-bfgs'.format(outputsuffix,phylipfile)
    if '-T' not in raxml_options:
        raxml_options += ' -T {0}'.format(multiprocessing.cpu_count())
    # Run raxml on renamed phylip
    run_raxml(phylipfile,raxml_options) 
    # Rename all raxml output files
    rename_list = glob('RAxML_*' + outputsuffix)
    for f in rename_list:
        rename_sequences(f, mapping)

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'inputfasta',
        help='Fasta file to run through raxml'
    )

    parser.add_argument(
        'raxmloptions',
        default='',
        help='What options to pass on to raxml command. Don\'t specify -n or -s or --no-bfgs'
    )

    parser.add_argument(
        '--output_name',
        dest='outputsuffix',
        default='raxmlrun',
        help='What to use for the -n option of raxml[Default: %(default)s]'
    )

    return parser.parse_args()

def run_raxml(phyml, raxmloptions):
    '''
    Run raxml on a phymlfile
    '''
    # Run raxml
    raxml_cmd = 'raxml {0}'.format(raxmloptions)
    print raxml_cmd
    try:
        subprocess.check_call(raxml_cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print '!!! There was an error running raxml !!!'
        print
        print e
        print
        print '----------------'
