from functools import partial
from itertools import dropwhile, takewhile, starmap
from toolz.itertoolz import  accumulate, first, last
from toolz.functoolz import complement, compose
import sh
import re

plus_n = lambda x, n : x + n
zero = lambda x, y : x
# changing plus_n to alter strings won't quite work for the string transformations because
# adding gaps is dependent on there being a difference between ref and pos.
dx_dy={
    'M' : (plus_n, plus_n), 'X' : (plus_n, plus_n), '=' : (plus_n, plus_n),
    'D' : (plus_n, zero), 'plus_n' : (plus_n, zero),
    'I' : (zero, plus_n), 'S' : (zero, plus_n),
    'H' : (zero, zero)
}

def transform(seqs, positions):
    (ref, query), ((rpos, qpos), (a, b)) = seqs, positions
    def fill(s, pos): return s[:pos] + '-'*abs(a-b) + s[pos:]
    if a < b: ref = fill(ref, rpos)
    elif a > b: query = fill(query, qpos)
    return  ref, query

def next_cigar_position(positions, cigarval):
    (count, sym), (refpos, qpos) = cigarval, positions
    f = dx_dy[sym]
    return f[0](refpos, count), f[1](qpos, count)
cig_diffs = partial(map, partial(next_cigar_position, (0, 0)))

def undo_cigar(x, y, z):
    cig = get_cig(z) #return reduce(transform, zip(gen_cigar(z), cig_diffs(z)), (x, y))
    return reduce(transform, zip(_gen_cigar(cig), cig_diffs(cig)), (x, y))

#in order to *totally* recreate the reads, you need to account for the other weird tag


split_cig = re.compile(r'(?:([0-9]+)([A-Z]))').findall
get_cig = lambda C: map(lambda x: (int(x[0]), x[1]), split_cig(C))
_gen_cigar = lambda c: accumulate(next_cigar_position, c, (0, 0))
gen_cigar = compose(_gen_cigar, get_cig)


def intersection(pred, seq): return last(takewhile(complement(pred), seq))
def firstwhere(pred, seq): return first(dropwhile(complement(pred), seq))
#reduce(next_cigar_position, cig_elements, (0,0))

# read = lambda x: vim.current.buffer.append(str(x))

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
['AAGC', 'AGTT', '1M1D1M1X1I']])) ==  [('AGG', 'A-G'), ('AGG', 'A-G'), ('T-T', 'TAT'), ('T---A', 'TATAA'), ('AAGC-', 'A-GTT')]
assert list(gen_cigar("1M2D1M")) == [(0, 0), (1, 1), (3, 1), (4, 2)]
assert intersection(lambda x: x[0] > 2, gen_cigar("1M2D1M")) == (1, 1)
assert firstwhere(lambda x: x[0] >= 2, gen_cigar("1M1D1M")) ==  (2, 1)
