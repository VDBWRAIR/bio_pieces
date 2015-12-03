'''
Command-line utility for querying VCF files. By default, outputs a full vcf file matching the query.

Usage:
    vcfcat filter <FILE1> ( --tag=<TAG> (--ge | --le | --gt | --lt | --eq | --ne) <VALUE> ) [-c]
    vcfcat exists <FILE1> (--tag=<TAG> ) [-c]
    vcfcat ambiguous <FILE1>  [-c]
    vcfcat vcall <FILE1>  [--csv | -c ]
    vcfcat diff <FILE1> <FILE2> (--tag=<TAG> [--threshold=<THRESH>])  [-c]
    vcfcat stat <FILE1>
    vcfcat statdiff <FILE1> <FILE2>

Options:
    -t=<THRESH> --threshold=<THRESH>  Numerical threshold to determine difference [Default: 0]
    -c, --count  Output the number of records instead of the actual file.
    --csv        Output variant calls in simplified tabular format
    --gt         Get records Greater Than <VALUE>
    --ge         Greater than or Equal
    --lt         Get records Less Thaan <VALUE>
    --leq        Get records Less than or equal to <VALUE>
    --eq         Get records Exactly Equal to <VALUE>
    --ne        Get records Not Equal to <VALUE>

Arguments:
    filter:  print only vcf records matching the filter as a new vcf file.
    exists:  only those records where the tag has a value.
    ambiguous:  only those records where the Called Base is ambiguous.
    vcall:  those records where the called base and reference base are different.
    diff:   those records in <FILE1> which were different than the corresponding records (at
    the same position) in <FILE2> according to some criteria. Providing a --tag argument will show
    show positions where the files have different entries in that field. The threshold argument (an integer)
    is used to determine if the values are "different enough", i.e., must be 10 more, or 10 less, etc.
    stat:   print simple information (like mean, etc.) for all VCF fields in a tabular format.
    statdiff:   print the (arithmetic) difference between the two file's "stat" output.
'''


from schema import Schema, Use,  Optional
from docopt import docopt
import operator
from bio_bits import vcfcat
import sys
import vcf
#TODO: Find a better way to dispatch commandline apps
ops = ['--ge', '--le', '--gt' , '--lt' , '--eq' , '--ne']
def validate_value(val):
    if val is None or not val.isdigit():
        return val
    else:
        return int(val)

def validate_vcf_file(name):
    if name is None:
        return None
    f  = open(name)
    return list(vcf.Reader(f))
#    return And( Use(open, error='could not open vcf file,'),
#               Use(vcf.Reader, error='pyvcf could not read file'),
#               Use(list))

def run(raw_args):
    '''
    Validates the arguments by converting the VCF files into lists of VCFRecords, converting
    threshold to int, and calling validate_value on <VALUE>.
    returns the command string for dispatching; the new args dictionary after validating, and the first
    filename so that that file can later be used for a vcf.VCFWriter object.
    :param dict raw_args: the result of docopt(__doc__)
    :return str, dict, str: the command string, the new args dict, and the first filename.
    '''
    commands = ['vcall', 'filter', 'diff', 'stat', 'statdiff', 'exists', 'ambiguous']
    schema_dict=    {
            '<FILE1>' : Use(validate_vcf_file),
            Optional('<FILE2>') : Use(validate_vcf_file),
            Optional('<VALUE>') : Use(validate_value),
            '--count' : bool,
            '--threshold' : Use(int, error='Threshold must be integer'),
            Optional('--tag') : lambda t: True,
            '--csv' : bool
             #tags.__contains__ #    error='Tag was not valid, should be one of {0}'.format(' '.join(tags)))
        }
    schema_dict.update( dict( (arg, bool) for arg in commands + ops))
    _args = Schema(schema_dict).validate(raw_args)
    cmd_str = [k for (k, arg) in _args.items() if k in commands and arg][0]
    filename = raw_args['<FILE1>']
    return cmd_str, _args, filename

def dispatch_cmd(cmd_str, args):
    '''
    Dispatch the appropriate vcfcat.py function matching the command string (the first commandline argument)
    :param str cmd_str: the command chosen, aka filter
    :param dict args: validated and type-converted command-line arguments
    :return list result: a list of VCFRecords matching the arguments
    '''
    cmd = getattr(vcfcat, cmd_str)
    if cmd_str == 'filter':
       op = [k for k, arg in args.items() if k in ops and arg ][0]
       opfunc = getattr(operator, op[2:])
       result = vcfcat.filter(args['<FILE1>'], args['--tag'], opfunc, args['<VALUE>'])
    elif cmd_str in ['stat', 'ambiguous']:
        result = cmd(args['<FILE1>'])
    elif cmd_str == 'diff':
        result = vcfcat.diff(args['<FILE1>'], args['<FILE2>'], args['--tag'], args['--threshold'])
    elif cmd_str == 'vcall':
        result = cmd(args['<FILE1>'])
    elif cmd_str == 'exists':
        result = cmd(args['<FILE1>'], args['--tag'])
    else: #'statdiff'
        result = vcfcat.statdiff(args['<FILE1>'], args['<FILE2>'])
    return result

def compute():
    raw_args = docopt(__doc__, version='Version 0')
    cmd_str, args, filename = run(raw_args)
    result = dispatch_cmd(cmd_str, args)
    return result, args, cmd_str, filename


def print_variant_call(rec):
    '''
    print vcalls in a simple tabular format
    '''
    print("{0}\t{1}\t{2}\t{3}".format(
        rec.CHROM, rec.POS, rec.REF, rec.INFO['CB']
    ))

def main():
    result, args, cmd, filename = compute()
    template = vcf.Reader(open(filename))
    stdout = sys.stdout
    if args['--count']:
        print(str(len(result)))
    elif args['--csv']:
        print('Reference\tPosition\tReference Base\tCalled Base')
        for rec in result:
            print_variant_call(rec)
    elif cmd.startswith('stat'):
        print(str(result))
    else:
        writer = vcf.Writer(stdout, template)
        for rec in result:
            writer.write_record(rec)

    #TODO: Make sure can accept input from stdin

if __name__ == "__main__":
    main()
