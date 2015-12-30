from __future__ import print_function
from os.path import dirname,basename,abspath,exists
import os
from functools import partial
import argparse
import shlex
import subprocess
import sys
from bio_bits.compat import OrderedDict, open

import sh

# Users cannot have these in the other args passed
STATIC_BLAST_ARGS = [
    '-num_threads', '-db', '-query'
]

# Users cannot have these in the other args passed
STATIC_DIAMOND_ARGS = [
    '-p', '--threads', '-d', '--db', '-q', '--query', '--daa', '-a'
]

# Static args for parallel
PARALLEL_ARGS = ['-u', '--pipe', '--cat', '--block', '10', '--recstart', '>', '--round-robin']

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
    if set(STATIC_BLAST_ARGS).intersection(shlex.split(blastoptions)):
        raise ValueError("You cannot supply any of the arguments inside of {0} as" \
            " optional arguments to blast".format(STATIC_BLAST_ARGS))
    args = list(PARALLEL_ARGS)
    args += generate_sshlogins(ninst)
    blast_path = sh.which(blasttype)
    blast_cmd = [blast_path]
    if blast_path is None:
        raise ValueError("{0} is not in your path(Maybe not installed?)".format(
            blasttype
        ))
    if task is not None:
        blast_cmd += ['-task', task]
    blast_cmd += ['-db', db,]
    blast_cmd += [blastoptions]
    blast_cmd += ['-query', '{}']
    args += [' '.join(blast_cmd)]
    cmd = sh.Command('parallel')
    run(cmd, *args, _in=open(inputfile), _out=open(outfile,'w'))

def parallel_diamond(inputfile, outfile, ninst, db, task, diamondoptions):
    '''
    Runs diamond commands in parallel on a given fasta file

    Will not run more than 1 diamond process per host as diamond utilizes
    threads better than blast

    Since diamond v0.7.9 produces a daa file, diamond view is required to output
    the tsv format that is similar to blast's output format. diamond view is
    automatically called on the produced .daa file so that GNU Parallel can combine
    all output into a single stream.

    :param str inputfile: Input fasta path
    :param str outfile: Output file path
    :param int ninst: number of threads to use if not in PBS or SGE job
    :param str db: Database path to blast against
    :param str task: blastx or blastp
    :param str diamondoptions: other options to pass to blast
    '''
    if set(STATIC_DIAMOND_ARGS).intersection(shlex.split(diamondoptions)):
        raise ValueError("You cannot supply any of the arguments inside of {0} as" \
            " optional arguments to diamond".format(STATIC_DIAMOND_ARGS))
    # This seems kinda stupid that we are just replacing cpu count for each
    # node with 1, but it is easier than refactoring other code to be better
    sshlogins = generate_sshlogins(ninst)
    for i in range(0,len(sshlogins),2):
        cpu,host = sshlogins[i+1].split('/')
        sshlogins[i+1] = '1/{0}'.format(host)
    dmnd_path = sh.which('diamond')
    if dmnd_path is None:
        raise ValueError("diamond is not in your path(Maybe not installed?)")
    # Diamond base command arguments
    # parallel replaces {} with the temporary file it is using
    #  and replaces {#} with the current file segment it is using
    # After diamond is finished, diamond view will be used to output the tsv
    # format of the file
    diamond_cmd = [
        dmnd_path, task, '--threads', str(ninst), '--db', db, '--query', '{}',
        '--daa', '{}.{#}', ';', dmnd_path, 'view', '--daa', '{}.{#}.daa'
    ]
    if len(sshlogins) > 2:
        args = list(PARALLEL_ARGS)
        args += sshlogins
        diamond_cmd_str = ' '.join(diamond_cmd) + diamondoptions
        args += [diamond_cmd_str]
        cmd = sh.Command('parallel')
        run(cmd, *args, _in=open(inputfile), _out=open(outfile,'w'))
    else:
        dcmd = sh.Command('diamond')
        args = [task]
        if diamondoptions:
            args += shlex.split(diamondoptions)
        p = run(
            dcmd, *args, threads=ninst, db=db, query=inputfile, daa=outfile
        )
        p = run(
            dcmd, 'view', daa=outfile+'.daa', _out=open(outfile,'w')
        )

def run(cmd, *args, **kwargs):
    '''
    Runs and prints what is being run to stdout
    '''
    kwargsignore = ['_in', '_out']
    kwargs_str = ' '.join(['--'+a+' '+str(v) for a,v in kwargs.items() 
        if a not in kwargsignore])
    args_str = ' '.join(args)
    print("[cmd] {0} {1} {2}".format(cmd._path, args_str, kwargs_str))
    try:
        p = cmd(*args, **kwargs)
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
