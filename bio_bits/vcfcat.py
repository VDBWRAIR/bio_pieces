from __future__ import print_function
'''
Usage:
    vcfcat filter <FILE> [ --tag  <TAG> (--ge | --le | --geq | --leq | --eq | --neq) <VALUE> ]... [-c]
    vcfcat filter <FILE> (--tag <TAG> --exists | --ambiguous)
    vcfcat diff <FILE1> <FILE2> (--tag <TAG> [--threshold <THRESH>])  [-c]
    vcfcat stat <FILE>
    vcfcat statdiff <FILE1> <FILE2>

Options:
    -t THRESHOLD --threshold=THRESHOLD  Numerical threshold to determine difference [Default: 0]
    -c, --count  Output the number of records.
    --no-header  Output without the VCF header info.
    --

    -o DIR --output=DIR  The output directory [default: vcf_compare_output]

'''
'''
vcf.infos => ordered dict of {'TAG' : vcf.parser.Info}
vcf will automatically handle the type=Interger, etc. field for us.
Note that ALT might be multiple alleles
'''
import pandas as pd
import operator
from functools import partial


def match_records(left, right, tag, threshold):
    '''
    :param VCFRecord left:
    :param VCFRecord right:
    :param str tag: the field to compare on
    :param int threshold: the amount difference necessary to return False
    :return bool: False if they are different according to threshold
    '''
    l, r = flatten_vcf(left), flatten_vcf(right)
    #assert tag in l and tag in r, "Tag not found in flattened record {0}".format(str(flatten_vcf))
    if operator.xor( (tag in l), (tag in r)):
        return False
    if tag not in l and tag not in r:
        return True
    if not threshold:
        return l[tag] == r[tag]
    else:
        return abs(l[tag] - r[tag]) < threshold

def vcall(records):
    '''
    :param list records: list of VCFRecord objects
    :return list: list of VCFRecords where the reference does not match the called base
    '''
    return [rec for rec in records if rec.REF != rec.INFO['CB']]

def diff(left, right, tag, threshold):
    '''
    :param list left: vcf Records
    :param list right: vcf Records
    :param str tag: i.e. CB, ALT
    :param int thershold: the amount of difference necessary
    :return list of tuples from left, right where the values differed
    '''
    #return [(l, r) for l, r in zip(left, right) if not match_records(l, r, tag, threshold)]
    return [l for l, r in zip(left, right) if not match_records(l, r, tag, threshold)]

def compare_value(comparator, val, rec, tag):
    '''
    :param builtin_function_or_method comparator: i.e. operator.eq, operator.lt
    :param object val: the object rec[tag] will be compared against
    :param VCFRecord rec:
    :param str tag: the tag we're interested in the rec object, i.e. 'CBD'
    :return bool: the result of the comparison
    '''
    flat_vcf = flatten_vcf(rec)
    assert tag in flat_vcf, "Tag not found in flattened record {0}".format(str(flatten_vcf))
    if type(flat_vcf[tag]) != type(val):
        raise ValueError("Type of record field {0}, {1}, does not match type of val {2}, {3}"
                         .format(flat_vcf[tag], type(flat_vcf[tag]), val, type(val)))
    return comparator(flat_vcf[tag], val)

def not_members(vcf_list, tag, collection):
    '''
    Get a list comprehension of all records in vcf_list where their field(tag) is not inside collection.
    :param list vcf_list: a list of VCFRecords
    :param str tag: the field of interest
    :param list collection: a list of objects
    :return list: all records where record[tag] or record.INFO[tag] is  not in collection
    '''
    return [rec for rec in vcf_list if (tag not in flatten_vcf(rec)) or flatten_vcf(rec)[tag] not in collection]

ambiguous = partial(not_members, tag='CB', collection=['A', 'C', 'G', 'T'])
exists = partial(not_members, collection= [None, [None], '-'])


def filter(vcf_list, tag, opfunc, value):
    '''
    i.e.  CBD > 12, etc.
    :param list vcf_list: list of vcf.model._Record objects as returned from vcf.Reader
    :param str operator: a function contained in the `operator` python stdlib module
    :param str tag: a valid field in the vcf file (i.e. ALT, CBD)
    :param object value: str, int/float, or list. the value to run operator against
    :return a subset of a vcf record where the operator evaluates to true
    '''
    #ALTs are stored as lists by vcf
    if tag == 'ALT' and type(value) is str:
        value = [base.strip()  for base in value.split(',')]
    if type(opfunc) is str:
        opfunc = getattr(operator, opfunc)
    check_value = partial(compare_value, opfunc, value)
    return [rec for rec in vcf_list if (tag in flatten_vcf(rec)) and check_value(rec, tag)]

def flatten_list(A):
  return A[0] if type(A) == list else A

def flatten_vcf(record):
    '''
    :param VCFRecord record: VCFRecord object
    :return dict: all fields as a flat dictionary, including those fields in rec.INFO
    '''
    _ids = ['CHROM', 'REF', 'QUAL', 'POS', 'ALT', 'FILTER', 'FORMAT', 'ID', 'INFO']
    fields = [getattr(record, _id) for _id in _ids]
    d = dict( (_id, value) for  _id, value in zip(_ids, fields) )
    d.update(dict((_id, flatten_list(field)) for _id, field in record.INFO.items()))
    return d

def vcf_file_to_df(vcf_records):
    '''
    Convert a list of vcf Records to a pandas DataFrame.
    '''
    return pd.DataFrame(flatten_vcf(rec) for rec in vcf_records)

def stat(records):
    '''
    Get simple statistics for each field like mean, total, count, etc.
    '''
    return vcf_file_to_df(records).describe()

def statdiff(recs_a, recs_b):
   '''
   Get the difference of the simple statistics using DataFrame.describe()
   '''

   df1, df2 = vcf_file_to_df(recs_a), vcf_file_to_df(recs_b)
   return df1.describe() - df2.describe()
