from os import environ
from sys import stderr
from time import time


def timeit(method):
    if "POLYA_BENCHMARK" not in environ:
        return method

    def timed(*args, **kw):
        ts = time()
        result = method(*args, **kw)
        te = time()

        stderr.write(f"[PERFORMANCE] {method.__name__} {int(te - ts)}s\n")

        return result

    return timed
