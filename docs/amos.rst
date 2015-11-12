====
AMOS
====

AMOS is a file format that is similar to any assembly file format such as ACE or SAM.
It contains information about each read that is used to assemble each contig.

The format is broken into different message blocks. For the Ray assembler, it 
produces an AMOS file that is broken into 3 types of message blocks

* RED

    .. code-block:: none

        {RED
        iid:\d+
        eid:\d+
        seq:
        [ATGC]+
        .
        qlt:
        [A-Z]+
        }

    iid
      Integer identifier
    eid
      Same as iid?
    seq
      Sequence data
    qlt
      Should be quality, but is only a series of D's from Ray assembler

* TLE

    .. code-block:: none

        {TLE
        src:\d+
        off:\d+
        clr:\d+,\d+
        }

    src
      RED iid that was used
    off
      One would think offset, but unsure what it actually means
    clr
      Not sure what this is either

* CTG

    .. code-block:: none

        {CTG
        iid:\d+
        eid:\w+
        com:
        .*$
        .
        seq:
        [ATGC]+
        .
        qlt:
        [A-Z]+
        .
        {TLE
        ...
        }
        }

    iid
      integer id of contig
    eid
      contig name
    com
      Communication software that generated this contig
    seq
      Contig sequence data
    qlt
      Supposed to be contig quality data, but for Ray it only produces D's
    TLE
      0 or more TLE blocks that represent RED sequences that compose the contig

Parsing
-------

bio_bits contains an interface to parse a given file handle that has been opened
on an AMOS file.

To read in the AMOS file you simply do the following

.. code-block:: python

    from bio_bits import amos
    a = None
    with open('AMOS.afg') as fh:
        a = amos.AMOS(fh)

CTG
^^^

To get information about the contigs(CTG) you can access the ``.ctgs`` attribute.
The contigs are indexed based on their iid so to get the sequence of contig iid 1 
you would do the following:

.. code-block:: python

    ctg = a.ctgs[1]
    seq = ctg.seq

To retrieve all the reads(RED) that belong to a specific contig:

.. code-block:: python

    reads = []
    for tle in ctg.tlelist:
        reads.append(a.reds[tle.src])

RED
^^^

To get information about the reads(RED) you can access the ``.reds`` attribute.
The reds are indexed based on their iid so to get the sequence of red iid 1 you
would do the following:

.. code-block:: python

    red = a.reds[1]
    seq = red.seq

If you want to convert a RED entry into anything you can use the ``.format``
method. The ``.format`` method allows you to utilize any of the properties of
a RED object such as ``.iid``, ``.eid``, ``.seq``, ``.qlt``. You can see in
the examples below how to do this.

Examples
--------

Here is an example of how to convert all RED blocks into a single fastq file

.. code-block:: python

    from bio_bits import amos

    # Fastq format string
    fastq_fmt = '@{iid}\n{seq}\n+\n{qlt}'

    with open('amos.fastq','w') as fh_out:
        with open('AMOS.afg') as fh_in:
            for iid, red in amos.AMOS(fh_in).reds.items():
                fq = red.format(fastq_fmt)
                fh_out.write(fq + '\n')
