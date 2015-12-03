#!/usr/bin/env python

# Designed to allow checkpointing of BEAST output files
# Input: beast_analysis.xml output_file_1.log output_file_2.log ... tree_file.trees
# Expect at least 3 files
# This matches column names in log files to parameter names in the original XML to set initial conditions

from __future__ import print_function

import sys
import re
import fileinput
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import argparse
import pandas as pd

def parse_logfile(logfilepath):
    '''
    Parse beast .log file and return pandas data frame

    :param str logfilepath: path to log file
    :param str indexcol: column to index(passed to pd.read_csv)
    :return: pandas DataFrame
    '''
    df = pd.read_csv(logfilepath, sep='\t', comment='#')
    return df

def hash_last_log_entries(*logfiles):
    '''
    Create a dictionary keyed by each unique log file header column name
    with the values from the last log entry 

    :param list logfiles: list of logfile paths
    :returns: dict of param values
    '''
    param = OrderedDict()
    # Iter each logfile
    for path in logfiles:
        # Load each file
        df = parse_logfile(path)
        # Last entry in current logfile
        lastentry = df.tail(1)
        _thing = None
        for i, row in lastentry.iterrows():
            _thing = OrderedDict(row)
        # Update param dict with dataframe as dictionary
        param.update(_thing)
    return param

def parse_args():
    parser = argparse.ArgumentParser('Restart beast from existing files')
    parser.add_argument(
        'xml_filename',
        help='Original xml input file path'
    )
    parser.add_argument(
        'tree_filename',
        help='Tree output file'
    )
    parser.add_argument(
        'log_filenames',
        nargs='+',
        help='Log files produced'
    )
    return parser.parse_args()

def main():
    #xml_filename = sys.argv[1]
    #tree_filename = sys.argv[-1]
    #log_filenames = sys.argv[2:-1]
    args = parse_args()
    xml_filename = args.xml_filename
    tree_filename = args.tree_filename
    log_filenames = args.log_filenames

    # Go through log files and create a hash of parameter name to value
    param = hash_last_log_entries(*log_filenames)

    # If there is a string of names in the form name1 \t name2 \t name3, this is a parameter array
    # and its map should be name to a list of values separated by spaces
    additional = OrderedDict()
    for name, value in param.items():
        m = re.findall('^(\S+[^0-9])([0-9]+)$', name)
        if m:
            shortname = m[0][0]
            #value = m[0][1]
            if shortname in additional:
                additional[shortname] += " {0}".format(value)
            else:
                additional[shortname] = "{0}".format(value)
    param.update(additional)

    sys.stderr.write('Param and values that will be used from log file\n')
    sys.stderr.write('{0}\n'.format(param))

    # Remove non-initialized parameters from param.
    if 'treeModel.rootHeight' in param:
        del param['treeModel.rootHeight']

    # Go through trees file and construct hash of id to taxa name
    taxa = OrderedDict()
    trees = open(tree_filename).readlines()
    for line in trees:
        line = line.strip()
        m = re.search("(\d+)\s+'*([A-Za-z0-9\-\_\.\/\| ]+)'*", line)
        if m:
            m = m.groups()
            taxa[m[0]] = m[1]
        if 'tree STATE_' in line:
            break


    tree = trees[-1]
    if not tree.startswith('tree STATE'):
        tree = trees[-2]
    # Go through last tree and replace ids with taxa names
    for taxid, taxa in taxa.items():
        tree = re.sub('(?<=[(,])'+taxid+'(?=\[|:)', "'"+taxa+"'", tree)
    # Remove tree STATE_... from tree
    tree = re.sub('tree STATE_.+ = \[\&R\] ', '', tree)
    # Remove [&rate=1.0] from tree(beast complained about it hopefully not needed)
    tree = tree.replace('[&rate=1.0]', '')

    # Go through XML file and replace initial parameter values with values from param
    # Only output lines if set to true
    outputline = True
    for line in fileinput.input(xml_filename):
        line = line.rstrip()
        # Don't output coalescentTree
        if '<coalescentTree id="startingTree"' in line:
            outputline = False

        m = re.search('parameter id="([^"]+)"', line)
        if m:
            name = m.group(1)
            #print 'Found line with param id now looking for param["{0}"]'.format(name)
            value = param.get(name)
            if value:
                #print 'Looking to replace {0} with {1} in this line'.format(name,value)
                if 'value=' in line:
                    line = re.sub('value="([^"]+)"', 'value="{0}"'.format(value), line)
                else:
                    line = re.sub('parameter id="([^"]+)"', 'parameter id="{0}" value="{1}"'.format(name,value), line)

        if outputline:
            print(line)

        if '</coalescentTree>' in line:
            outputline = True
            print('\t<newick id="startingTree">')
            print('\t\t'+tree)
            print('\t</newick>')

if __name__ == '__main__':
    main()
