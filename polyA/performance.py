from os import environ
from time import time


def timeit(method):
    if "POLYA_BENCHMARK" not in environ:
        return method

    def timed(*args, **kw):
        ts = time()
        result = method(*args, **kw)
        te = time()

        print(f"[BENCHMARK] {method.__name__} {int(te - ts)}s")

        return result

    return timed
