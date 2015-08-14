aln_pca
=======

aln_pca is used to build a pca 3D plot for a given fasta alignment file.
The image produced will contain 9 3D plots of the same plot just rotated slightly
so you can visualize each axis easier.

By default aln_pca will build the simplest identity matrix by comparing each
sequence and assigning the value 1 to each matching nucleotide and 0 to any 
mismatch and then summing the values to produce the pairwise identity for each
sequence.

Running aln_pca on test dataset
-------------------------------

Default substitution matrix
+++++++++++++++++++++++++++

Here we are running the simple pairwise identity matrix(default) to generate
our output file at docs/_static/pca.png

.. code-block:: bash

    $> aln_pca tests/testinput/aln1.fasta --outfile docs/_static/pca.png

Produces:

.. image:: /_static/pca.png

Providing substitution matrix file
++++++++++++++++++++++++++++++++++

Now we will provide the included jalview_snm.txt substitution matrix to generate
the output file at docs/_static/jal_view.png

.. code-block:: bash

    $> aln_pca tests/testinput/aln1.fasta --substitution-matrix jalview_snm.txt --outfile docs/_static/jalview.png

Produces:

.. image:: /_static/jalview.png

Substitution Matrix
-------------------

Here we will use a very simple aligned fasta file

Input fasta alignment::
    
    >id1  
    ATGC
    >id2
    AAGC
    >id3
    CGTR
    >id4
    ACGT

This would produce the following identity matrix by default:

+---+---+---+---+---+
|   |id1|id2|id3|id4|
+===+===+===+===+===+
|id1|  4|  3| 0 |  2|
+---+---+---+---+---+
|id2|  3|  4|  0|  2|
+---+---+---+---+---+
|id3|  0|  0|  4|  0|
+---+---+---+---+---+
|id4|  2|  2|  0|  4|
+---+---+---+---+---+


By default aln_pca will use the simplest substitution matrix where for every 
non-substitution a 1 is used and any substitution yields a 1

If you specify a substitution matrix file on the command line you can
change how the identity matrix is generated.

We can use the SNS substitution matrix used by `Jalview`_::

         A   C   G   I   N   R   T   U   X   Y  ~
    A   10  -8  -8  1   1   1   -8  -8  1   -8  1
    C   -8  10  -8  1   1   -8  -8  -8  1   1   1
    G   -8  -8  10  1   1   1   -8  -8  1   -8  1
    I   1   1   1   10  1   0   1   1   0   0   1
    N   1   1   1   1   10  1   1   1   1   1   1
    R   1   -8  1   0   1   10  -8  -8  0   -8  1
    T   -8  -8  -8  1   1   -8  10  10  1   1   1
    U   -8  -8  -8  1   1   -8  10  10  1   1   1
    X   1   1   1   0   1   0   1   1   10  0   1
    Y   -8  1   -8  0   1   -8  1   1   0   10  1
    ~   1   1   1   1   1   1   1   1   1   1   1

*Note*: Any non listed base will be assigned a 1 as it will be looked up using the
``~`` key

This matrix is also included with the source code of this project as 
``jalview_snm.txt``

Since this is just a space separated file we can just paste it into a file and
supply the file path using the ``--substitution-matrix`` option as you can see
above in the examples.

Any tab delimited file will work as long as it has the Nucleotides for the X and Y
axis and has numerical values for each cell. If your matrix does not include a
lookup for every substitution possible it will error or you can include the 'all'
lookup ``~`` as the table above does. What this does is for any base that is not
listed in the table, it will instead lookup using the ``~`` listed value. For
example, if we use the above substitution table and tried looking up a substitution
of an ``A`` to a ``W`` it would yield the value ``1``.

If we use the same input fasta sequence as above using the JalView substitution
matrix, it would yield the following identity matrix:

+---+---+---+---+---+
|   |id1|id2|id3|id4|
+===+===+===+===+===+
|id1| 40| 22|-32|  4|
+---+---+---+---+---+
|id2| 22| 40|-32|  4|
+---+---+---+---+---+
|id3|-32|-32| 40|-32|
+---+---+---+---+---+
|id4|  4|  4|-32| 40|
+---+---+---+---+---+

Output file
-----------

Once the program completes you will will find a file called ``pca.png`` that
contains the 9 plots stacked vertically.

You can optionally specify where to put the output file using the ``--outfile``
argument to ``aln_pca``


Why 9 plots?
++++++++++++

Since this is a 3D plot, it was decided to view the plot from 9 different viewpoints.
Essentially the plot is viewed with no tilt at 0, 45 and 90 degrees rotated.

Then the plot is tilted 45 degrees down and rotated again at 0, 45 and 90 degrees.

The last time the plot is tilted 90 degrees(looking straight down) and rotated
at 0, 45 and 90 degrees.

.. _jalview: http://www.jalview.org/help/html/calculations/scorematrices.html#simplenucleotide
