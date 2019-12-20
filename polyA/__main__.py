import sys

from ._cli import Options, run

if __name__ == "__main__":
    options = Options(sys.argv[1:])
    run(options)