from logging import Logger
from os import environ
from sys import stderr
from time import time
from typing import Optional, TextIO

ENV_VAR_NAME = "POLYA_PERFORMANCE"


def timeit(method, file: TextIO = stderr, logger: Optional[Logger] = None):
    """
    A benchmarking helper intended for us as a decorator. To use,
    decorate a function with `@timeit`, then set the `POLYA_PERFORMANCE`
    environment variable to any value when the program is run.

    If a logger is provided to the decorator then it will be used,
    otherwise output will be written to stderr unless a `file` is provided,
    in which case that will be used.

    TODO: Allow time units to be specified at decoration site

    >>> import os, io
    >>> stream = io.StringIO()
    >>> os.environ[ENV_VAR_NAME] = 'true'
    >>> perf_sum = timeit(sum, file=stream)
    >>> perf_sum([1, 2])
    3
    >>> stream.getvalue()
    '[PERFORMANCE] sum 0s\\n'
    """
    if ENV_VAR_NAME not in environ:
        return method

    def timed(*args, **kw):
        ts = time()
        result = method(*args, **kw)
        te = time()

        msg = f"[PERFORMANCE] {method.__name__} {int(te - ts)}s\n"

        if logger is None:
            file.write(msg)
        else:
            logger.info(msg)

        return result

    return timed
