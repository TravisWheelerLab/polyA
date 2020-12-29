# A quick script to simplify profiling with Scalene.
# Use this script as the Python file to be executed
# (profiled).
#
# $ scalene <scalene options> profile.py <polya options>

from polyA._app import run
run()

