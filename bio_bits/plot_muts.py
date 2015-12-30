'''
Usage:
    plot_muts.py --query <query> --refs <refs> [--out <outfile>]

Options:
    --refs,-r=<refs>     Fasta file, sequence with earliest year is base reference
    --query,-q=<query>   Query sequences
    --out,-o=<outfile>   Figure saved here

Help:
    All sequences must be the same length.
'''
from __future__ import print_function
import matplotlib.patches as mpatches
import numpy as np
import range_regex
from functools import partial
import operator
import os, sys, re
from Bio import SeqIO
import matplotlib.pyplot as plt
import docopt, schema
from operator import itemgetter as get
try:
    #below import is necessary for some reason
    from scipy.stats import poisson
    import scipy
    DISTRIBUTION = scipy.stats.poisson
except ImportError:
    DISTRIBUTION = None
compose2 = lambda f,g: lambda *x: f(g(*x)) 
compose = lambda *fs: reduce(compose2, fs)

years = range_regex.range_regex.regex_for_range(1900, 2015)
year_regex = re.compile(years)
hamming = compose(sum, partial(map, operator.ne))
legend = {"queries": 'r', "references": 'b', "interval": 'g'}
'''it seems like pdist gives results that are too small to be useful?'''
#def pdist(s1, s2):
#    assert len(s1) == len(s2), "All sequences must be the same length! %s %s" % (s1, s2)
#    return hamming(s1, s2)/float(len(s1))

def extract_year(header):
    #s = header[-4:]
    s = header.split('/')[3]
    try: 
        return int(year_regex.search(s).group())
    except:
        print( header)
# had to add 2015 to A/England/50220895/


def get_seqs_and_years(fn):
    fasta = SeqIO.parse(fn, format="fasta")
    info = [ (str(seq.seq), seq.id) for seq in fasta]
    seqs, ids = zip(*info)
    years = map(extract_year, ids)
    return seqs, years


def process(refs_fn, query_fn, save_path=None):
    ref_seqs, ref_years = zip(*sorted(zip(*get_seqs_and_years(refs_fn)), key=get(1)))
    ref_seqs = map(str.upper, ref_seqs)
    super_ref_seq, super_ref_year = ref_seqs[0], ref_years[0]
    get_mutations = partial(hamming, super_ref_seq)
    def get_relative_info(seqs, years):
         muts = map(get_mutations, seqs)
         dists = [yr - super_ref_year for yr in years]
         return muts, dists
    ref_muts, ref_dists =  get_relative_info(ref_seqs[1:], ref_years[1:])
    query_muts, query_dists = get_relative_info(*get_seqs_and_years(query_fn))
    do_plot(ref_dists, ref_muts, query_dists, query_muts, save_path)
    #map(compose(print, '{0}\t{1}'.format ), ref_dists, ref_muts)

def do_plot(x1, y1, x2, y2, save_path=None):

    ax = plt.subplot(111)
    max_x = max(max(x1), max(x2))
    plot_muts(ax, x1, y1, color=legend['references'], polyfit=False, max_x=max_x)
    plot_muts(ax, x2, y2, color=legend['queries'])
    legend_info = [mpatches.Patch(label=n, color=c) for n, c in legend.items()]
    """ http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot"""
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(handles=legend_info, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel("Years since Base reference")
    plt.ylabel("p-distance")
    if save_path:
        plt.savefig(save_path)
    plt.show()

def plot_muts(ax, x, y, color, dist=DISTRIBUTION, polyfit=False, max_x=None):
    #problem was didn't account for +b
    '''if norm distribution, probably have to scale (via passing loc= and scale=)'''
    ax.scatter(x, y, color=color)
    if polyfit:
        ''' this forces a polyfit with y-intercept at zero, necessary because
        we necessarily start with 0 mutations from the query at year 0.'''
        x = np.array(x)[:,np.newaxis]
        m, _, _, _ = np.linalg.lstsq(x, y)
        x, y = np.linspace(0,max_x, 100), m*np.linspace(0,max_x, 100)
        ax.plot(x, y, color='y')
    if dist:
        """see  http://stackoverflow.com/a/14814711/3757222"""
        R = dist.interval(0.95, y)
        interval_left, interval_right = R
        interval_color = legend['interval']
        ax.plot(x, interval_left, color=interval_color)
        ax.plot(x, interval_right,color=interval_color)

#def test_more():
#    refs = range(25), range(25)
#    queries = [1, 5, 20, 10], [2, 20, 40, 10]
#    do_plot(refs[0], refs[1], queries[0], queries[1], None)
#
#def test_plot():
#    ''' can verify this works by using scipy.stats.norm.interval instead'''
#    default_x = range(25)
#    default_y = range(0, 50, 2)
#    plot_muts(default_x, default_y, 'r', True, scipy.stats.poisson, max_x=max(default_x))
#    plt.show()

def main():
    #if sys.argv[1] == 'test': test_more()
    scheme = schema.Schema(
        { '--query' : os.path.isfile,
          '--refs' : os.path.isfile,
         schema.Optional('--out') : lambda x: True
        # schema.Or(lambda x: x is None,  #check file can be created
        #                                      lambda x: os.access(os.path.dirname(x), os.W_OK))
         })
    args = docopt.docopt(__doc__, version='Version 1.0')
    scheme.validate(args)
    queries, refs, out = args['--query'], args['--refs'], args['--out']
    process(refs, queries, out)

if __name__ == '__main__': main()

'''
'Our Data Example'!A1:A2991,'Our Data Example'!B1:B2991, 'Our Data Example'!D1:D2991, 'Our Data Example'!E1:E2991
'Our Data Example'!B1:B8, 'Our Data Example'!D1:D8, 'Our Data Example'!E1:E8

'Our Data Example'!A1:A999,'Our Data Example'!B1:B999, 'Our Data Example'!D1:D999, 'Our Data Example'!E1:E999

'''
