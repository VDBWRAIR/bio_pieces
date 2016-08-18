'''
Usage:
  genbank <infile>
'''
from Bio import SeqIO
import dateparser
import sys
from docopt import docopt

def rec_qualifier(rec, key):
  for f in rec.features:
    if f.qualifiers:
      if 'collection_date' in f.qualifiers:
        return f.qualifiers[key]

def fix_date(d):
    date = dateparser.parse(d)
    return date.strftime('%Y-%m-%d')

def rec_to_fasta(rec):
  raw_collection_date = rec_qualifier(rec, 'collection_date')[0]
  country = rec_qualifier(rec, 'country')[0]
  date = dateparser.parse(raw_collection_date)
  date = date.strftime('%Y-%m-%d')
  rec.description = "|{}|{}".format(country, date)
  return rec

  # acc = rec.name
  #header = "{}|{}|{}".format(acc, country, date)
  #seq = str(rec.seq)
  #return rec.format("fasta")
def run(infile):
  with open(infile) as input:
    gb = SeqIO.parse(input, 'genbank')
    results = map(rec_to_fasta, gb)
    SeqIO.write(results, sys.stdout, 'fasta')

def main():
  raw_args = docopt(__doc__, version='Version 1.0')
  run(raw_args['<infile>'])
  sys.exit(0)


# ['17-Feb-2010']
# rec = degen.id_to_record("KJ627355")
