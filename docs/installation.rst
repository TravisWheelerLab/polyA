Installation
============
PolyA may be consumed as an ordinary Python package:

Install using the `pip` tool:

    ``pip install polyA``

Docker
--------------
There is a `Dockerfile` in the repo root. The image it describes is
used for running tests in CI and can be used locally for convenience.

When dependencies change the image must be rebuilt and the new
version pushed to Docker Hub. This can be done with `make container`
if you have the correct permissions. Otherwise, ask a maintainer
to do it for you.

There is also a file called `Dockerfile_run`. This is a runner image that is
available for users who prefer to run PolyA through Docker (which is often more
convenient). The runner image can be built with `tool/build-runner-image.sh`
and, assuming sufficient permissions, pushed to Docker Hub with
`tool/push-runner-image.sh`.

The scripts assume you are building the latest version of the image and so set
the `latest` tag. If this is not the case, it is fairly simple to issue the
correct Docker commands manually.

A runner image can be found on `<https://hub.docker.com/repository/docker/traviswheelerlab/polya>`_


External software dependencies
=================================

Easel
--------------
PolyA internally adjusts alignment scores to account for the scale multiplier
implicit to each score matrix (the so-called `lambda` value).
Unless the lambda value is given in the matrix file, PolyA must compute it; this
is done using the `esl_scorematrix` utility found in the
[Easel](https://github.com/EddyRivasLab/easel/) software package.

To prepare that tool, follow these steps:
::

    # get source code
    % cd $HOME/git  # or wherever you like to place repositories
    % git clone https://github.com/EddyRivasLab/easel
    % cd easel
    % autoconf
    % ./configure
    % make
    % gcc -g -Wall -I. -L. -o esl_scorematrix -DeslSCOREMATRIX_EXAMPLE esl_scorematrix.c -leasel -lm


    # add the easel directory to your path, e.g.
    % export PATH=$HOME/git/easel:$PATH
    # you may wish to add this to your path permanently, e.g.
    #    https://opensource.com/article/17/6/set-path-linux

    # alternatively, you can set the --easel-path argument to $HOME/git/easel

ULTRA
----------
PolyA can optionally include tandem repeat annotations from ULTRA as competitors in
the adjudication process. Options are:

(i) If you have an ultra annotation output on hand for the sequence being annotated, reference
the ultra output file with the --ultra-data flag, e.g.

``% polyA --ultra-data ultra.file ALIGNMENTS MATRICES``

(ii) If you wish to run [ULTRA](https://github.com/TravisWheelerLab/ultra) on the fly, you'll
need the software:
::

    % cd $HOME/git  # or wherever you like to place repositories
    % git clone https://github.com/TravisWheelerLab/ULTRA ultra
    % cd ultra
    % cmake .
    % make

then you'll reference the ultra binary in the PolyA command:

``% polyA --ultra-path $HOME/git/ultra/ultra --sequences seq.fasta  ALIGNMENTS MATRICES``


Development
=================
This project uses [Pipenv](https://pipenv.pypa.io/en/latest/), which can be
installed through Homebrew for Mac users. It must be installed before the
Makefile targets, or the other commands listed in this document will work.

In order to run a command which relies on the project virtual environment, such
as `python foo.py`, it is necessary to either run ``pipenv shell`` first, which
will put you into a shell that has the correct Python in its `PATH`, or prefix
the command with ``pipenv run`` (e.g. ``pipenv run python foo.py``).

To build a PyPI package, use ``make build-package``. To publish, use
``make publish-package``. We use [Flit](https://flit.readthedocs.io/en/latest/) to
handle building and publishing a package, so it is also possible to invoke Flit
directly to do anything else it is capable of.

Makefile
--------------
A Makefile is available, run ``make`` or ``make help`` in the project root to see
the available targets.

Dependencies
------------------
If you prefer not to use the Makefile, or if you need to add or remove
dependencies, the following commands will allow you to manage dependencies.

::

    # Get to work (from within project directory)
    # This drops you into a shell inside the project virtual environment which means
    # that commands that "pipenv run" may be elided from other commands.
    pipenv shell

    # Fetch runtime dependencies
    pipenv install

    # Fetch runtime and development dependencies
    pipenv install --dev

    # Add a runtime dependency
    pipenv install <package>

    # Add a development dependency
    pipenv install --dev <package>