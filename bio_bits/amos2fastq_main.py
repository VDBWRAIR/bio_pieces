'''
Usage: amos2fastq <fastqs>... --amos=<amos>

Given any number of fastq files and a specified AMOS file (usually with .afg extension) describing a series of contigs, create seperate fastq files arranged by the contig each read aligned to.

'''
from schema import Schema, Use, And
from docopt import docopt
from bio_bits import amos2fastq
#Do file validation immediately when script is started
def all_elemnts_unique(collection):
    return len(collection) == len(set(collection))

# possibly try/catch
def validate_args(raw_args):
    scheme = Schema({'<fastqs>' : And( [Use(open, error='fastq file not readable')],
                                 all_elemnts_unique, error='fastq files must be different and exist'),
        '--amos' : Use(open, error='AMOS file should be readable') })
    return scheme.validate(raw_args)

def main():
    raw_args = docopt(__doc__, version='Version 0')
    parsed_args = validate_args(raw_args)
    return amos2fastq.make_fastqs_by_contigs(parsed_args['<fastqs>'], parsed_args['--amos'])


if __name__ == "__main__":
    main()


