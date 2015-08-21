from __future__ import print_function
from os.path import dirname,basename,abspath,exists
import os
from functools import partial
import argparse
import shlex
import subprocess
import sys
try:
    from __builtin__ import open
except ImportError:
    from builtins import open

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import sh

# Staticly set options for blast
MAX_TARGET_SEQS = 10
BLAST_FORMAT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'inputfasta',
        help='fasta sequence file to blast as -query'
    )
    parser.add_argument(
        'outfile',
        help='blast output file to write that will contain all results'
    )
    parser.add_argument(
        '--ninst',
        type=int,
        default=1,
        help='Number of total cpus to use[Default: %(default)s]'
    )
    parser.add_argument(
        '--db',
        help='Blast db path'
    )
    parser.add_argument(
        '--blast_exe',
        choices=('diamond','blastn','blastx'),
        help='Executable to call'
    )
    parser.add_argument(
        '--task',
        choices=('megablast','dc-megablast','blastn','blastx','blastp'),
        default=None,
        help='task to use for blast/diamond'
    )
    parser.add_argument(
        '--blast_options',
        default='',
        help='Options to pass on to blast/diamond[Default: %(default)s]'
    )
    return parser.parse_args()

def parallel_blast(inputfile, outfile, ninst, db, blasttype, task, blastoptions):
    '''
    Runs blast commands in parallel on a given fasta file

    :param str inputfile: Input fasta path
    :param str outfile: Output file path
    :param int ninst: number of cpus to use if not in PBS or SGE job
    :param str db: Database path to blast against
    :param str blasttype: Blast exe to use(blastn, blastx, blastp)
    :param str task: Blast task to run with -task option for blasttype or 
        None if blastx/blastp
    :param str blastoptions: other options to pass to blast
    '''
    blast_path = sh.which(blasttype)
    if blast_path is None:
        raise ValueError("{0} is not in your path(Maybe not installed?)".format(
            blasttype
        ))
    args = ['-u', '--pipe', '--block', '10', '--recstart', '>']
    args += generate_sshlogins(ninst)
    args += [blast_path]
    if task is not None:
        args += ['-task', task]
    args += ['-db', db, '-max_target_seqs', str(MAX_TARGET_SEQS),
        '-outfmt', '"'+BLAST_FORMAT+'"'
    ]
    args += shlex.split(blastoptions)
    args += ['-query', '-']
    cmd = sh.Command('parallel')
    run(cmd, args, inputfile, outfile)

def parallel_diamond(inputfile, outfile, ninst, db, task, diamondoptions):
    '''
    Runs diamond commands in parallel on a given fasta file

    Will not run more than 1 diamond process per host as diamond utilizes
    threads better than blast

    :param str inputfile: Input fasta path
    :param str outfile: Output file path
    :param int ninst: number of threads to use if not in PBS or SGE job
    :param str db: Database path to blast against
    :param str task: blastx or blastp
    :param str diamondoptions: other options to pass to blast
    '''
    '''
    diamond -task blastx -compress 0 -db /path/nt -o outfile -query inputfile -o outfile
    my $cmd = "$type $task_option  $options -q  $query -d $db  -o $out"; 
    '''
    # This seems kinda stupid that we are just replacing cpu count for each
    # node with 1, but it is easier than refactoring other code to be better
    sshlogins = generate_sshlogins(ninst)
    for i in range(0,len(sshlogins),2):
        cpu,host = sshlogins[i+1].split('/')
        sshlogins[i+1] = '1/{0}'.format(host)
    dmnd_path = sh.which('diamond')
    if dmnd_path is None:
        raise ValueError("diamond is not in your path(Maybe not installed?)")
    args = ['-u', '--pipe', '--block', '10', '--recstart', '>', '--cat']
    args += sshlogins
    args += [
        dmnd_path, task, '--threads', str(ninst), '--db', db, '--query', '{}',
        '--compress', '0', '-a', outfile
    ] + shlex.split(diamondoptions)
    cmd = sh.Command('parallel')
    run(cmd, args, inputfile, outfile)

def run(cmd, args, infile, outfile):
    print("[cmd] {0} {1}".format(cmd._path, ' '.join(args)))
    try:
        p = cmd(*args, _in=open(infile), _out=open(outfile,'w'))
        print(p)
    except sh.ErrorReturnCode as e:
        print("There was an error")
        print(str(e.stderr))
        sys.exit(e.exit_code)

def get_hostfile():
    '''
    Return the path to the hostfile|machinefile depending on SGE or PBS
    If local, returns ''
    '''
    return os.environ.get('PBS_NODEFILE', os.environ.get('PE_HOSTFILE', ''))

def parse_hostfile(hostfile_fh):
    '''
    Parse job queue system's hostfile for either PBS style or SGE style

    These are available with $PBS_NODEFILE and $PE_HOSTFILE during jobs that are run
    Always assumes first column is hostname
    If second column exists, then it is cpu count
    Any other columns are ignored

    :param file hostfile_fh: file like object that contains host file information
    '''
    hosts = OrderedDict()
    for line in hostfile_fh:
        x = line.split()
        hostname = x[0]
        if hostname not in hosts:
            hosts[hostname] = 0
        cpus = 1 # Default is to increment 1 for ever found(machinefile)
        if len(x) > 1: #Non PBS machinefile
            try:
                cpus = int(x[1])
            except ValueError as e:
                raise ValueError(
                    "Invalid second column value found in hostfile: {0}".format(x[1])
                )
        hosts[hostname] += cpus
    return list(hosts.items())

def generate_sshlogins(ninst=None):
    '''
    Generate a list of --sshlogins compatible with GNU Parallel such that they
    match PBS_NODEFILE, PE_HOSTFILE or just the local host

    :param int ninst: Number of cpus to use if not in SGE or PBS job
    :returns: ['--sshlogin cpu/host1', ..., '--sshlogin cpu/hostN']
    '''
    path = get_hostfile()
    sshlogins = []
    if not path:
        if ninst is None:
            ninst = ''
        else:
            ninst = '{0}/'.format(ninst)
        sshlogins = ['--sshlogin', '{0}:'.format(ninst)]
    else:
        with open(path) as fh:
            hosts = parse_hostfile(fh)
        if len(hosts) == 1: # Single Node
            sshlogins = ['--sshlogin', '{0}/:'.format(hosts[0][1])]
        else: # Multiple nodes
            for x in hosts:
                sshlogins.append('--sshlogin')
                sshlogins.append('{1}/{0}'.format(*x))
    return sshlogins

def main():
    args = parse_args()
    assert exists(args.inputfasta), '[error] {0} does not exist'.format(args.inputfasta)
    if args.blast_exe == 'diamond':
        parallel_diamond(
            args.inputfasta, args.outfile, args.ninst, args.db, args.task,
            args.blast_options 
        )
    else:
        parallel_blast(
            args.inputfasta, args.outfile, args.ninst, args.db,
            args.blast_exe, args.task, args.blast_options 
        )
