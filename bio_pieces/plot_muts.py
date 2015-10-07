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
import matplotlib.patches as mpatches
import numpy as np
import range_regex
from funcy import compose
from functools import partial
import operator
import os, sys, re
from Bio import SeqIO
import matplotlib.pyplot as plt
import docopt, schema
from operator import itemgetter as get
#below import is necessary for some reason
from scipy.stats import poisson
import scipy
import itertools
#from scipy.stats import norm

years = range_regex.range_regex.regex_for_range(1900, 2015)
year_regex = re.compile(years)
hamming = compose(sum, partial(map, operator.ne))
def pdist(s1, s2):
    assert len(s1) == len(s2), "All sequences must be the same length! %s %s" % (s1, s2)
    return hamming(s1, s2)/float(len(s1))
#pdist is hamming divided by lenght of sequence
def extract_year(s): return int(year_regex.search(s).group())
def get_seqs_and_years(fn):
    fasta = SeqIO.parse(fn, format="fasta")
    info = [ (str(seq.seq), seq.id) for seq in fasta]
    seqs, ids = zip(*info)
    years = map(extract_year, ids)
    return seqs, years

legend = {"queries": 'r', "references": 'b', "interval": 'g'}

def process(refs_fn, query_fn, save_path=None):
    ref_seqs, ref_years = zip(*sorted(zip(*get_seqs_and_years(refs_fn)), key=get(1)))
    super_ref_seq, super_ref_year = ref_seqs[0], ref_years[0]
    get_mutations = partial(pdist, super_ref_seq)
    def get_relative_info(seqs, years):
         muts = map(get_mutations, seqs)
         dists = [yr - super_ref_year for yr in years]
         return muts, dists
    ref_muts, ref_dists =  get_relative_info(ref_seqs[1:], ref_years[1:])
    query_muts, query_dists = get_relative_info(*get_seqs_and_years(query_fn))
    do_plot(ref_dists, ref_muts, query_dists, query_muts, save_path)

def do_plot(x1, y1, x2, y2, save_path):
    max_x = max(max(x1), max(x2))
    plot_muts(x1, y1, color=legend['references'], interval=True, polyfit=True, max_x=max_x)
    plot_muts(x2, y2, color=legend['queries'], interval=False)
    legend_info = [mpatches.Patch(label=n, color=c) for n, c in legend.items()]
    plt.legend(handles=legend_info)
    plt.xlabel("Years since Base reference")
    plt.ylabel("p-distance")
    if save_path:
        plt.savefig(save_path)
    plt.show()

def plot_muts(x, y, color, interval=False, dist=scipy.stats.poisson, polyfit=True, max_x=None):
    plt.scatter(x, y, color=color)
    if polyfit:
        m = np.polyfit(x, y, 1)[0]
        x, y = np.linspace(0,max_x), m*np.linspace(0,max_x)
        plt.plot(x, y, color='y')
    if interval:
        ''' can verify this works by using scipy.stats.norm.interval instead.'''
        """see  http://stackoverflow.com/a/14814711/3757222"""
        R = dist.interval(0.95, y)
        interval_left, interval_right = R
        interval_color = legend['interval']
        plt.plot(x, interval_left, color=interval_color)
        plt.plot(x, interval_right,color=interval_color)

def test_plot():
    ''' can verify this works by using scipy.stats.norm.interval instead'''
    default_x = range(25)
    default_y = range(0, 50, 2)
    plot_muts(default_x, default_y, 'r', True, scipy.stats.poisson, max_x=max(default_x))
    plt.show()

def main():
    if sys.argv[1] == 'test': test_plot()
    scheme = schema.Schema(
        { '--query' : os.path.isfile,
          '--refs' : os.path.isfile,
         schema.Optional('--out') : lambda x: True
        # schema.Or(lambda x: x is None,  #check file can be created
        #                                      lambda x: os.access(os.path.dirname(x), os.W_OK))
         })
    args = docopt.docopt(__doc__, version='Version 1.0')
    #args = scheme.validate(raw_args)
    scheme.validate(args)
    queries, refs, out = args['--query'], args['--refs'], args['--out']
    process(refs, queries, out)

if __name__ == '__main__': main()
