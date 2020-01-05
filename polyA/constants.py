import math
import sys
import numpy as np

"""
The width of the "window" used to break a sequence up into pieces.
Measured in base pairs.
"""
CHUNK_SIZE = 30

"""
The string used to represent "not a number" for serialization and
deserialization.
"""
NAN_STRING = "NaN"

"""
The base probability that we do not change our estimate from one
base pair to the next. This will generally be a very small value.
"""
SAME_PROB = float(np.nextafter(1, 0))

# +-------------------------------------------------------------+
# | Computed constants, do not manually modify below this line. |
# +-------------------------------------------------------------+

"""
The base probability of changing our estimate from one family of
sequences to another, part way through adjudicating a sequence.
"""
CHANGE_PROB = 1.0 - SAME_PROB

LOGGED_CHANGE_PROB = math.log(CHANGE_PROB)
# LOGGED_CHANGE_PROB = -106.10123583452  # Perl value override

LOGGED_SAME_PROB = math.log(SAME_PROB)
# LOGGED_SAME_PROB = 0  # Perl value override
