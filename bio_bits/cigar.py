from functools import partial
from itertools import dropwhile, takewhile, starmap, combinations, imap, ifilter
import itertools
from toolz.itertoolz import  accumulate, first, last, drop, mapcat, take
from toolz.functoolz import complement, compose
from operator import itemgetter, add
import fn
import sh
import re 
import multiprocessing  as mp
import sys

#TODO: in order to *totally* recreate the reads, you need to account for the MD (mismatching positions) tag.  
#This requires running bwa with different options. http://bio-bwa.sourceforge.net/bwa.shtml

#def ilen(seq): sum(1 for _ in seq)
#indel_pos = compose(ilen, partial(takewhile, '-'.__ne__))

# adding gaps is dependent on there being a difference between ref and pos.

def transform(seqs, positions):
    ''' given ref & query and positions and offsets, builds a sort of alignment by filling in 
    gaps appropriately with '-'.
    :param tuple seqs reference, query strings
    :param tuple positions  
    :return tuple of reference, query strings ''' 
    (ref, query), ((rpos, qpos), (a, b)) = seqs, positions
    def fill(s, pos): return s[:pos] + '-'*abs(a-b) + s[pos:]
    if a < b: ref = fill(ref, rpos)
    elif a > b: query = fill(query, qpos)
    return  ref, query

no_add = lambda x, y : x
dx_dy={
    'M' : (add, add), 'X' : (add, add), '=' : (add, add),
    'D' : (add, no_add), 'add' : (add, no_add),
    'I' : (no_add, add), 'S' : (no_add, add),
    'H' : (no_add, no_add)
}
def next_cigar_position(positions, cigarval, withtype=False):
    '''compute the next position for the reference and position by "applying"
    the position change that letter indicates, using `dx_dy` defined above.
    if withtype is True, also returns the letter itself, and expects the last count
    and last letter..'''
    if withtype: (count, sym), (refpos, qpos, oldcount, oldsym) = cigarval, positions
    else: (count, sym), (refpos, qpos) = cigarval, positions
    f = dx_dy[sym]
    if withtype: return f[0](refpos, count), f[1](qpos, count), count, sym
    return f[0](refpos, count), f[1](qpos, count)

split_cig = re.compile(r'(?:([0-9]+)([A-Z]))').findall
cig_diffs = partial(map, partial(next_cigar_position, (0, 0)))
get_cig = lambda C: map(lambda x: (int(x[0]), x[1]), split_cig(C))
_gen_cigar = lambda c: accumulate(next_cigar_position, c, (0, 0))
gen_cigar = compose(_gen_cigar, get_cig)

def undo_cigar(x, y, z):
    cig = get_cig(z) 
    return reduce(transform, zip(_gen_cigar(cig), cig_diffs(cig)), (x, y))

def get_insert_quals(ref, pos, bases):

    '''separate i/o and destructuring from logic.
    make generic for indel/snp.'''
    view = sh.samtools('view', "{ref}:{pos}-{pos}")
    #split tabs, etc.
    has_insert = lambda l: 'I' in l[5]
    # etc.
    #samtools calcmd?

assert undo_cigar('AAGC', 'AGTT', '1M1D1M1X1I') == ('AAGC-', 'A-GTT')
assert list(starmap(undo_cigar, [
["AGG", "AG"   ,"1M1D1M"],
["AGG", "AG"   ,"1M1D1M"],
["TT" , "TAT"  , "1M1I1M"],
["TA" , "TATAA", "1M3I1M"],
['AAGC', 'AGTT', '1M1D1M1X1I']])) ==  \
        [('AGG', 'A-G'), ('AGG', 'A-G'), ('T-T', 'TAT'), ('T---A', 'TATAA'), ('AAGC-', 'A-GTT')]
assert list(gen_cigar("1M2D1M")) == [(0, 0), (1, 1), (3, 1), (4, 2)] 

def intersection(pred, seq): return last(takewhile(complement(pred), seq))
def firstwhere(pred, seq): return first(dropwhile(complement(pred), seq))

assert intersection(lambda x: x[0] > 2, gen_cigar("1M2D1M")) == (1, 1)
assert firstwhere(lambda x: x[0] >= 2, gen_cigar("1M1D1M")) ==  (2, 1)


assert reduce(next_cigar_position,  \  #corresponds to ('AAGC-', 'A-GTT')
       takewhile(lambda x: x[1] != 'I', get_cig('1M1D1M1X1I')), (0, 0)) == (4, 3)

def mutpos(cig_str, types):
    ''':return (refpos, mutpos'''
    assert any(map(types.__contains__, cig_str)), "Cigar string %s does not have mutation of type %s" % (cig_str, types)
    return reduce(next_cigar_position, takewhile(lambda x: x[1] not in types, get_cig(cig_str)), (0, 0))

indelpos = partial(mutpos, types='DI')
query_indel_pos = compose(last, indelpos)
ref_indel_pos = compose(first,  indelpos)

ref_del_pos = compose(first, partial(mutpos,  types='D'))
ref_ins_pos = compose(first, partial(mutpos,  types='I'))
ref_snp_pos = compose(first, partial(mutpos,  types='X'))

assert ref_ins_pos('1M1D1M1X1I') == 4

def mutations(cig_str):
    return drop(1, accumulate(partial(next_cigar_position, withtype=True), get_cig(cig_str), (0, 0, 0, 0)))

#TODO: Unfortnuately, need to use the MD tag, because bwa does not provide
# the 'X' tag for mismatches.
snps = compose(partial(filter, lambda x: x[-1] == 'X'), mutations)

def snp_rpos_base(qstring, cig_str):
    return compose(partial(mapcat, partial(rpos_base, qstring)), snps)(cig_str)

def rpos_base(qstring, (rpos, qpos, count, type)):
    ''':return refpos, snp'''
    return ((qstring[(qpos-count)+i], rpos +i) for i in xrange(0, count)) 

try:
    mutpos('1M1D1M', 'X')
except AssertionError:
    pass

def get_info(samrow):
    #TODO: extract query_string and cigar_string from file
    qname, mpos, cigar, qstring = itemgetter(0, 3, 5, 9)(samrow.split('\t'))
    mpos = int(mpos)
    _snps = snp_rpos_base(qstring, cigar)
    #have to add map positions to computed rpos's.
    return map(lambda (b,i):(b, i+mpos), _snps)

''' :lhn int minimum shared mutations in order to mark a group of reads as a local haplotype cluster.
    :ghn int minimum shared mutations in order to mark a group of reads as a GLOBAL haplotype cluster.  '''
if __name__ == '__main__':
    '''split by RG (read group), then get_info, want to keep track of qname to prob.
       Find all local and global haplotype clusters.'''
    lhn = ghn = 3
    shareCount = compose(len, lambda (x,y): x & y) 
    is_hap = lambda n: lambda x: shareCount(x) > n
    _fn = sys.argv[1]
    p = mp.Pool()
    shortReads = imap(get_info, take(3, sh.samtools.view(_fn, _iter=True)))
    shortReads = imap(set, shortReads)
    cartProduct = itertools.product(*list(itertools.tee(shortReads)))
    shortReads = ifilter(is_hap(lhn), cartProduct)
    # shortReads  should become the cluster of local haplotypes.
    for r in shortReads:
        if r: print r 

#ghs = filter(is_hap(ghn), combinations(longReads, longReads))
#lhs = filter(is_hap(lhn), combinations(shortReads, shortReads)) 
