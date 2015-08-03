#!/usr/bin/env python
# encoding: utf-8
"""
ctleptop.py

Created by Dereje Jima on May 21, 2015
"""
import time
import sys
import os.path
from Bio.Seq import *
from Bio.Alphabet import IUPAC
from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
from itertools import groupby
from Bio.Data import CodonTable
from Bio.Data.IUPACData import ambiguous_dna_values
import yaml
import argparse
__docformat__ = "restructuredtext en"

AMBICODON = {"R": ["A", "G"], "Y": ["C", "T"], "W": ["A", "T"], "S": ["G", "C"], "K": ["T", "G"],
                                                "M": ["C", "A"], "D": ["A", "T", "G"], "V": ["A", "C", "G"], "H": ["A", "C", "T"],
                                                "B": ["C", "G", "T"], "N" : ["A", "C", "T", "G"]
                                                }
#def readFasta(fasta_name):
    #"""string -> tuple
    #Given a fasta file. Yield tuples of header, sequence
    #:fasta_name: Name of a file to read
    #:returns: tuple of fasta header and sequence line

    #"""
    #fh = open(fasta_name)
    ## ditch the boolean (x[0]) and just keep the header of sequence since we
    ## know they alternate
    #fasta_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    #for header in fasta_iter:
        ## drop the ">"
        #header = header.next()[1:].strip()
        ## join all sequence line to one
        #seq = "".join(s.strip() for s in fasta_iter.next())
        #yield header, seq


def getNearbyChars(nt): return AMBICODON.get(nt) or nt

def nearbyPermutations(letters,index=0):
    if (index >= len(letters)): return set([''])
    subWords = nearbyPermutations(letters, index + 1)
    nearbyLetters = getNearbyChars(letters[index])
    return permutations(subWords, nearbyLetters)

def  permutations(subWords, nearbyLetters):
   permutations = set()
   for subWord in subWords:
       for letter in nearbyLetters:
           permutations.add(letter + subWord)
   return permutations


def getaalist(codonlist):
    """Convert codon list to aa.

    :codonlist: a list of codon
    :aalist: return aa list

    """
    aalist = []
    for codon in codonlist:
            aa = Seq(codon, IUPAC.unambiguous_dna)
            aa = str(translate(aa))
            aalist.append(aa)
    return aalist


def list_overlap(list1, list2):
    """(str, list) -> bool
    Return True  if the two list hava element that overlaps.

    >>>list_overlap('RAC',['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y'])
    True
    >>>list_overlap('ACT',['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y'])
    False

    """
    for i in list1:
        if i in list2:
            return True
    return False

# define our method
def replace_all(text, dic):
        for i, j in dic.iteritems():
            text = text.replace(i, j)
        return text


def access_mixed_aa(file_name):
    """(str) ->(list,list,list,list).
    Return four lists by reading multifasta nucleotide fasta file
    >>>access_mixed_aa(file_name)

    """
    from Bio import SeqIO
    aa = []
    nucleotide_idx=[]
    nucl_codon = []
    seqids=[]
    for seq_record in SeqIO.parse(file_name, 'fasta'):
        seq_id = seq_record.id
        seq_len= len(seq_record)
        header,seqline=seq_record.id, str(seq_record.seq)
    #for header, seqline in readFasta(file_name):
        # print header + "\n" + seq_line

        #my_seq = Seq(seq_line, IUPAC.extended_dna)
        my_seq = Seq(str(seqline), IUPAC.ambiguous_dna)
        #seq2 = Seq("ARAWTAGKAMTA", IUPAC.ambiguous_dna)
        #seq2 = seq2.translate()
        # print seq2
        # print ambiguous_dna_values["W"]
        # print IUPAC.ambiguous_dna.letters
        seqline = seqline.replace("-", "N")
        n = 3
        codon_list = {i+n: seqline[i:i + n]  for i in range(0, len(seqline), n)}
        #print yaml.dump(ambi_codon)
        #print yaml.dump(codon_list)
        ambi_nucl = AMBICODON.keys()
        #print ambi_nucl
        # print ambi_codon["Y"]
        for key, codon in sorted(codon_list.iteritems()):
            #print "key: ", key , "codon:", codon
            if list_overlap(codon, ambi_nucl):
                d, e, f = codon
                m = [d, e, f]
                #print codon, ".....", key
                # print type(ambi_nucl)
                items = set(m).intersection(ambi_nucl)
                indexm = m.index(list(items)[0])
                #print "index ...", indexm
                items = list(items) # eg. ['R']
                for idx, val in enumerate(items):
                    # print idx
                    # print val
                    item = items[idx]

                    #print item
                    #print ambi_codon[item]
                    codonlist=list(nearbyPermutations(codon))
                    #print "codon list :", codonlist
                    val = getaalist(codonlist)
                    val = list(set(val)) # remove if aa codon is the same eg. ['D', 'D']
                    val = "/".join(val) # yeild 'I/L'
                    val = str(val)
                    #print "codonlist *****", codonlist
                    #print "aa val *******", val
                    if "/" in val and indexm == 2:
                        key = key
                        nucleotide_idx.append(key)
                        nucl_codon.append(codon)
                        seqids.append(seq_id)
                    elif "/" in val and indexm == 1:
                        key = key-1
                        nucleotide_idx.append(key)
                        nucl_codon.append(codon)
                        seqids.append(seq_id)
                    elif "/" in val and indexm == 0:
                        key = key-2
                        nucleotide_idx.append(key)
                        nucl_codon.append(codon)
                        seqids.append(seq_id)
                    else:
                        pass
                    #print ".....", val
                    aa.append(val)

            else:
                #print "codon3 ..." ,codon
                aa1 = Seq(codon, IUPAC.unambiguous_dna)
                aa1 = aa1.translate()
                aa1 = str(aa1)
                aa.append(aa1)
    return aa, nucleotide_idx, nucl_codon,seqids

def main():
    parser = argparse.ArgumentParser(description='Convert inframe nucleotide fasta file'\
                                     'to protein codes and report mixed aa with its locations')
    parser.add_argument("-i", type=str, help="input nucleotide fasta file")
    args = parser.parse_args()
    file_name = args.i
    #print list(nearbyPermutations("AAR"))
    #print listReplace("YAR", ["Y", "R"], {"Y": ['C', 'T'],  "R":["A", "G"]})
    my_list = access_mixed_aa(file_name)
    #print my_list
    aa,nuc_idx,nucl_codon,seqids = access_mixed_aa(file_name)
    #print aa,nuc_idx,nucl_codon,seq_id

    #print ''.join(aa)
    #print nuc_idx
    #print listReplace('GAY', 'Y', ['C', 'T'])
    #print listReplace("ARA", "R", ["A", "G"])
    #print listReplace("WTA", "W", ["A", "T"])

    import re
    #print aa[331]
    #print aa
    pattern = re.compile(r'.+\/.+')
    amb_aa_codon = []
    amb_aa_indx =[]
    for indx, letter in enumerate(aa):
        #print indx, ".....", letter
        if pattern.match(letter):
            amb_aa_codon.append(letter)
            amb_aa_indx.append(indx +1)
       #     print indx + 1, letter
    #print(amb_aa_codon)
    my_list = zip(seqids,nuc_idx,amb_aa_indx,nucl_codon,amb_aa_codon)
    #print my_list
    my_list = [list(elem) for elem in my_list]
    #print list(my_list)

    from tabulate import tabulate
    #print my_list
    print tabulate(my_list, headers=['seq id','nt Position','aa position', 'nt composition', 'aa composition'])



if __name__ == '__main__':
        main()
