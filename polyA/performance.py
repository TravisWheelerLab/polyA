from logging import Logger
from os import environ
from sys import stderr
from time import time
from typing import Optional


def timeit(method, logger: Optional[Logger] = None):
    if "POLYA_PERFORMANCE" not in environ:
        return method

    def timed(*args, **kw):
        ts = time()
        result = method(*args, **kw)
        te = time()

        msg = f"[PERFORMANCE] {method.__name__} {int(te - ts)}s\n"

        if logger is None:
            stderr.write(msg)
        else:
            logger.info(msg)

        return result

    return timed
