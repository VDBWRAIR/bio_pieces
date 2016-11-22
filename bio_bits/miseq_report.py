#!/usr/bin/env python

import gzip
import multiprocessing
import sys
import logging
import itertools

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from path import Path
import sh
import pandas as pd

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('miseqreport')

try:
    from itertools import imap
except ImportError:
    pass

class MiSeqRun(object):
    def __init__(self, runpath=None):
        self.runpath = Path(runpath)
        self.basecalls = self.runpath / 'Data' / 'Intensities' / 'BaseCalls'
        self.samplesheet = self.basecalls / 'SampleSheet.csv'
        self.rundate = self.runpath.basename().split('_')[0]

    def get_fastq_gz(self):
        return self.basecalls.glob('*_R1_*.fastq.gz') + \
            self.basecalls.glob('*_R2_*.fastq.gz')

    def get_paired_fastq_gz(self):
        gzs = zip(
            self.basecalls.glob('*_R1_*.fastq.gz'),
            self.basecalls.glob('*_R2_*.fastq.gz')
        )
        return gzs

    def get_samplenames(self):
        samplenames = []
        for line in self.samplesheet.lines(retain=False):
            if line.endswith(',,'):
                samplenames.append(line.split(',')[0])
        return samplenames

    def sample_from_fq(self, fqpath):
        return fqpath.basename().split('_')[0]

    def get_fqs_per_samplename(self):
        persample = []
        for pair in self.get_paired_fastq_gz():
            persample.append((self.sample_from_fq(pair[0]), pair))
        return persample

    def all_sample_read_stats(self, threads=None):
        '''
        Return [(samplename, sample_read_stats(samplenamepair)), ...] which contains stats for each sample
        '''
        pool = multiprocessing.Pool(processes=threads)
        return sorted(pool.map(stat, self.get_fqs_per_samplename()))

    def run_stats(self, allstats=None, threads=None):
        '''
        CSV output for the run
        :param iterable allstats: (samplename, (sample_read_stats(forward), sampleread_stats(reverse))
        '''
        if allstats is None:
            allstats = self.all_sample_read_stats(threads=threads)
        headers = "SampleName,TotalReads,FReads,FQual,FLen,FBases,RReads,RQual,RLen,RBases"
        rowfmt = ','.join(['{{{0}}}'.format(i) for i in range(len(headers.split(',')))])
        report = Path('report.{0}.csv'.format(self.rundate))
        report.write_text(headers + '\n')
        for sn, pair in allstats:
            f,r = pair
            row = (sn,f[0]+r[0]) + f + r
            report.write_text(rowfmt.format(*row) + '\n', append=True)
        return report

def run_fastqc(fqpaths, outdir=None, threads=None, tmpdir=None):
    '''
    Run fastqc on fqpath
    '''
    # Remove old one
    op = Path(outdir)
    op.rmtree_p()
    op.mkdir()
    print sh.fastqc(fqpaths, threads=threads, outdir=outdir, dir=tmpdir)
    # build fastqc html index
    index = op / 'index.html'
    html = ['<table>']
    for link in sorted(op.glob('*.html')):
        html.append('<tr>')
        html.append('<td><a href="{0}">{0}</a></td>'.format(link.basename()))
        html.append('</tr>')
    html.append('</table>')
    index.write_lines(html)

def stat(x):
    return (x[0], sample_read_stats(x[1]))

def sample_read_stats(paired):
    '''
    Return (read_stats(forward), read_stats(reverse))
    '''
    return map(read_stats, paired)

def read_stats(readpath):
    '''
    Return (#reads, avgqual, avglen, totalbases) for a fastq.gz read
    '''
    reads = 0
    quals, lengths = [],[]
    with gzip.open(readpath) as fh:
        logger.info('Getting stats for {0}'.format(readpath.basename()))
        for title, seq, qual in FastqGeneralIterator(fh):
            reads += 1
            lengths.append(len(seq))
            qual = map(lambda q: ord(q)-33, qual)
            quals.append(sum(qual)/float(len(qual)))
    totalbases = sum(lengths)
    return (reads, sum(quals)/float(len(quals)), totalbases/float(len(lengths)), totalbases)

def html_report(reportpath, fastqcdir):
    '''
    reportpath - path to report.csv
    fastqcdir - Path.py.Path object representing directory of fastqc outputs
    '''
    html_report_pth = reportpath + '.html'
    fastqc_htmls = sorted(map(lambda x: (x.basename().split('_')[0],x), fastqcdir.glob('*fastqc.html')))
    paired_htmls = dict(map(lambda x: (x[0], map(lambda j: j[1], x[1])), itertools.groupby(fastqc_htmls, lambda x: x[0])))
    # Colors
    defaultcolor = '#FFFFFF'
    # Color map for stddev +/- mean
    # 1 is the color for 1 stdev above mean
    # -1 is color for for 1 stdev below mean
    colors = {
        1: '#82E0AA',
        2: '#2ECC71',
        3: '#239B56',
        -1: '#C0392B',
        -2: '#922B21',
        -3: '#641E16'
    }
    # Read the csv report
    df = pd.read_csv(reportpath, index_col=0) 
    # Build table
    html = ['<table border="1">']
    # Put headers in table
    html.append('<tr>')
    html.append('<td>{0}</td>'.format(df.index.name))
    html += map(lambda h: '<td>{0}</td>'.format(h), df.keys())
    html.append('</tr>')
    # Put values in table
    for sn, values in df.iterrows():
        #print list(paired_htmls[sn])
        htmls = list(paired_htmls[sn])
        html.append('<tr>')
        # Samplename col
        href = '<a href="{0}">{1}</a>'
        fhref = href.format(htmls[0], 'R1')
        rhref = href.format(htmls[1], 'R2')
        html.append('<td>{0} ({1},{2})</td>'.format(sn,fhref,rhref))
        for col, val in values.iteritems():
            # The following two calculations should be moved out such that
            # they are only calculated once per column and then put into a dict
            # that can be looked up by colheader
            # Current column standard dev
            colstd = df[col].std()
            # Current column mean
            colmean = df[col].mean()
            # Creates a new mapping from colors based on stddev/mean
            cmap = map(lambda x: (colstd*x[0]+colmean, x[1]), colors.iteritems())
            # Find values above mean where val is greater than stddev
            colorv = filter(lambda x: x[0] >= colmean and val >= x[0], cmap)
            # If non found, then try looking below mean
            if colorv:
                color = max(colorv)[1]
            else:
                color = filter(lambda x: x[0] <= colmean and val <= x[0], cmap)
                if color:
                    color = min(color)[1]
                else:
                    color = defaultcolor
            html.append('<td bgcolor="{0}">{1}</td>'.format(color, val))
        html.append('</tr>')
    html.append('</table>')
    html_report_pth.write_lines(html)
    return html_report_pth

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print "Usage: report.py <miseqrunpath> <fastqcoutdir> <threads> <tmpdir>"
        sys.exit(1)
    runpath, fqcoutdir, threads, tmpdir = sys.argv[1:]
    threads = int(threads)
    run = MiSeqRun(runpath)
    logger.info("Creating csv report")
    reportpth = run.run_stats(threads=threads)
    logger.info("Wrote {0}".format(reportpth))
    logger.info("Creating {0}".format(fqcoutdir))
    logger.debug(str(run.get_fastq_gz()))
    run_fastqc(run.get_fastq_gz(), outdir=fqcoutdir, threads=threads, tmpdir=tmpdir)
    logger.info("Created {0}".format(fqcoutdir))
    logger.info("Writing html report")
    reportpth = html_report(reportpth, Path(fqcoutdir))
    logger.info("Wrote {0}".format(reportpth))
