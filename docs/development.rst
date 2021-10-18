Development
===========

This project uses `Pipenv <https://pipenv.pypa.io/en/latest>`_,
which can be installed through Homebrew for Mac users. It must
be installed before the Makefile targets, or the other commands
listed in this document will work.

In order to run a command which relies on the project virtual
environment, such as ``python foo.py``, it is necessary to
either run ``pipenv shell`` first, which will put you into a
shell that has the correct Python installed, or prefix the
command with ``pipenv run`` (e.g. ``pipenv run python foo.py``).

::

    # Get to work (from within project directory)
    # This drops you into a shell inside the project virtual
    # environment which means that commands that "pipenv run"
    # may be elided from other commands.
    pipenv shell

Makefile
--------

A Makefile is available, run ``make`` or ``make help`` in the
project root to see the available targets.

Dependencies
------------

The Makefile has targets for managing dependencies, but it may
be easier, particularly for Python developers, to simply invoke
Pipenv directly. Example commands are listed below.

If you prefer not to use the Makefile, or if you need to add
or remove dependencies, the following commands will allow you
to manage dependencies.

::

    # Fetch runtime dependencies
    pipenv sync

    # Fetch runtime and development dependencies
    pipenv sync --dev

    # Add a runtime dependency
    pipenv install <package>

    # Add a development dependency
    pipenv install --dev <package>

Unit Tests
----------

Unit tests use `pytest <https://pytest.org/en/latest/>`_. The
linked documentation contains examples of how to use it to
both write and run tests. The
`examples <https://pytest.org/en/latest/example/index.html>`_ are
quite extensive.

Documentation tests are also fine for simple cases.

::

    # Run tests
    pipenv run python -m pytest

    # or
    make check-fast check-slow

    # or
    make check

Docker Container
----------------

There is a ``Dockerfile`` in the repo root. The image it
describes is used for running tests in CI and can be used
locally for convenience.

When dependencies change the image must be rebuilt and the
new version pushed to Docker Hub. This can be done with
``make container`` if you have the correct permissions.
Otherwise, ask a maintainer to do it for you.

There is also a file called ``Dockerfile_run``. This is a
runner image that is available for users who prefer to run
PolyA through Docker (which is often more convenient). The
runner image can be built with ``tool/build-runner-image.sh``
and, assuming sufficient permissions, pushed to Docker Hub with
``tool/push-runner-image.sh``.

The scripts assume you are building the latest version of
the image and so set the ``latest`` tag. If this is not the
case, it is fairly simple to issue the correct Docker commands
manually.

This image can be found on Docker Hub:
https://hub.docker.com/repository/docker/traviswheelerlab/polya

Documentation
-------------

We use Sphinx for project documentation. Run ``make docs``
to update the documentation sources and build the HTML. Use
``make docs-serve`` to serve the documentation locally.

Release Process
---------------

Generating a new release is a four step process. The first step
MUST be completed first, but the rest may be done in any order.

First, update the value of ``VERSION`` in ``polyA/_version.py``,
incrementing the various portions of the version number
accordingly. We would like to follow
`Semantic Versioning <https://semver.org>`_, at least in general.
Commit your changes to the `master` branch.

Second, `create <https://github.com/TravisWheelerLab/polyA/releases/new>`_
a "release" on GitHub. The tag should consist of a "v", followed
immediately by the version number you chose above. For example:
``v1.2.3``.

Third, run ``make build-package``, followed by
``make publish-package`` assuming there aren't any build errors.

Fourth, and finally, create and push a new runner image to Docker
Hub. Build the image by running ``./tool/build-runner-image.sh``,
and publish it with ``./tool/push-runner-image.sh``. The version
tag will be set automatically and the ``latest`` tag will also be
updated.
