from logging import Logger
from os import environ
from sys import stderr
from time import time
from typing import Optional, TextIO

ENV_VAR_NAME = "POLYA_PERFORMANCE"

_logger: Optional[Logger] = None
_output_file: Optional[TextIO] = stderr


def set_timeit_file(file: TextIO):
    global _output_file
    _output_file = file


def reset_timeit_file():
    global _output_file
    _output_file = stderr


def set_timeit_logger(logger: Logger):
    global _logger
    _logger = logger


def reset_timeit_logger():
    global _logger
    _logger = None


def timeit(name: Optional[str] = None):
    """
    A benchmarking helper intended for us as a decorator. To use, decorate a
    function with `@timeit()`, then set the `POLYA_PERFORMANCE` environment
    variable to any value when the program is run.

    Output will be written to stderr unless a `file` is provided by calling
    `set_timeit_file`, in which case that will be used. If a logger is passed to
    `set_timeit_logger` then that will be used for output preferentially.

    >>> import os, io
    >>> stream = io.StringIO()
    >>> os.environ[ENV_VAR_NAME] = 'true'
    >>> set_timeit_file(stream)
    >>> perf_sum = timeit()(sum)
    >>> perf_sum([1, 2])
    3
    >>> stream.getvalue()[:31]
    '[PERFORMANCE] time spent in sum'
    """

    def decorator(method):
        if ENV_VAR_NAME not in environ:
            return method

        def timed(*args, **kw):
            ts = time()
            result = method(*args, **kw)
            te = time()

            if name is None:
                method_name = method.__name__
            else:
                method_name = name

            msg = f"time spent in {method_name} {float(te - ts):.4f}s"

            if _logger:
                _logger.info(msg)
            elif _output_file:
                _output_file.write(f"[PERFORMANCE] {msg}\n")

            return result

        return timed

    return decorator
