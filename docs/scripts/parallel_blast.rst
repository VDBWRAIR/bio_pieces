==============
parallel_blast
==============

Parallel blast is a wrapper script around the blast commands as well as diamond.
It utilizes GNU Parallel to run the commands in parallel by splitting up the input
fasta files and distributes them across multiple subprocesses. If it detects that
it is running inside of a PBS or SGE job it will run the job on multiple hosts
that may be allocated to the job.

parallel_blast requires that you have gnu parallel installed and in your environments
PATH as well as diamond and/or blastn/blastx/blastp.

* `diamond`_
* `blast`_
* `GNU parallel`_

Usage
=====

You can get all the arguments that can be supplied via the following

.. code-block:: bash

    $> parallel_blast --help

Examples
--------

For the examles below assume you have an input fasta in the current directory
called ``input.fasta``

Running blastn
++++++++++++++

.. code-block:: bash

    $> parallel_blast input.fasta output.blast --ninst 4 --db /path/to/nt \
    --blast_exe blastn --task megablast --blast_options "--evalue 0.01"
    [cmd] /path/to/parallel -u --pipe --block 10 --recstart > --sshlogin 4/: /path/to/blastn -task megablast -db /path/to/nt -max_target_seqs 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -query -

Notice how we had to quote the additional ``--blast_options``

Running diamond
+++++++++++++++

Diamond v0.7.9 is the version that was tested with parallel_blast. As diamond is
still in development the options may change and thus parallel_blast may not run
it correctly.

.. code-block:: bash

    $> parallel_blast input.fasta out.blast --ninst 4 --db /path/to/diamondnr \
    --blast_exe diamond --task blastx --blast_options "--tmpdir dtmp"
    [cmd] /path/to/parallel -u --pipe --block 10 --recstart > --cat --sshlogin 1/: /path/to/diamond blastx --threads 4 --db /path/to/diamondnr --query {} --compress 0 -a out.blast

Notice how even though we specified ``--ninst 4`` that ``--sshlogin 1/:`` was used
and ``--threads 4`` was set instead.

Command that is run
+++++++++++++++++++

You will notice in the examples above that when you run parallel_blast that it
outputs the command that it is running in case you want to copy/paste it and run
it yourself sometime.

You might notice that the command does not include all the quoted arguments such 
as the ``--recstart`` argument which should be ``--recstart ">"`` as well as 
the ``--outfmt`` which should be quoted as ``--outfmt "6 ..."``. If you intend on 
rerunning the command you will have to add the quotes manually.

Running inside of a PBS or SGE Job
==================================

parallel_blast is able to detect if it is running inside of a PBS or SGE job by
looking to see if ``PBS_NODEFILE`` or ``PE_HOSTFILE`` is set in the environment's
variables.

If it finds either of them it will run the job by supplying ``--sshlogin`` for each
host it finds in the file.

``PBS_NODEFILE`` and ``PE_HOSTFILE`` have different syntax so parallel_blast first
builds a CPU,NODENAME list from them.

PBS_NODEFILE
------------

This file is parsed and counts how many of each unique host is listed such that
the following PBS_NODEFILE::

    node1.localhost
    node2.localhost
    node2.localhost
    node3.localhost
    node3.localhost
    node3.localhost

would run 1 instance on node1.localhost, 2 instances on node2.localhost and 3
instances on node3.localhost

PE_HOSTFILE
-----------

This file is almost in the exact syntax that parallel_blast uses so it is almost
a 1-to-1 mapping.

Diamond and multiple hosts
--------------------------

Since diamond utilizes threads much more efficiently than blast, for each unique
host in a job only 1 instance is launched but the ``-p`` option is set to the number
of CPUS for each host listed in the ``PE_HOSTFILE`` or ``PBS_NODEFILE``

.. _diamond: https://github.com/bbuchfink/diamond
.. _blast: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST
.. _GNU parallel: http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
