Installation
============

PyPI
----

PolyA may be consumed as an ordinary Python package and
installed using the Pip tool:

::

    pip install polyA

    # Outside of a Virtual Environment, install as a user
    # package, which doesn't require root.
    pip install -U polyA

Docker
------

It is possible to run PolyA as a Docker container. This is
simpler to install as it includes all required Python code
and runtime dependencies.

The container image is on
`Docker Hub <https://hub.docker.com/repository/docker/traviswheelerlab/polya>`_.
See the README there for instructions on running it.

Runtime Dependencies
--------------------

Easel
^^^^^

PolyA internally adjusts alignment scores to account for the
scale multiplier implicit to each score matrix (the so-called
`lambda`` value). Unless the lambda value is given in the matrix
file, PolyA must compute it; this is done using the
``esl_scorematrix`` utility found in the
`Easel <https://github.com/EddyRivasLab/easel>`_ software
package.

To prepare that tool, follow these steps:

::

    # get source code
    cd $HOME/git  # or wherever you like to place repositories
    git clone https://github.com/EddyRivasLab/easel
    cd easel
    autoconf
    ./configure
    make
    gcc -g -Wall -I. -L. \
        -o esl_scorematrix \
        -DeslSCOREMATRIX_EXAMPLE \
        esl_scorematrix.c \
        -leasel -lm

    # add the easel directory to your path, e.g.
    export PATH=$HOME/git/easel:$PATH
    # you may wish to add this to your path permanently, e.g.
    # https://opensource.com/article/17/6/set-path-linux

    # alternatively, you can set the --easel-path argument to $HOME/git/easel

ULTRA
^^^^^

PolyA can optionally include tandem repeat annotations
from ULTRA as competitors in the adjudication process.
Options are:

(i) If you have an ultra annotation output on hand for
the sequence being annotated, reference the ultra output
file with the --ultra-data flag. For example:

::

    polyA --ultra-data ultra.file ALIGNMENTS MATRICES

(ii) If you wish to run
`ULTRA <https://github.com/TravisWheelerLab/ultra>`_ on
the fly, you'll need the software:

::

    cd $HOME/git  # or wherever you like to place repositories
    git clone git@github.com:TravisWheelerLab/ULTRA ultra
    cd ultra
    cmake .
    make

After that, you'll reference the ultra binary in the PolyA
command:

::

    polyA --ultra-path $HOME/git/ultra/ultra \
          --sequences seq.fasta  ALIGNMENTS MATRICES
