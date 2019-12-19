import math
import numpy as np

"""
The width of the "window" used to break a sequence up into pieces.
"""
CHUNK_SIZE = 30

"""
The base probability that we do not change our estimate from one base pair to
the next.
"""
SAME_PROB = float(np.nextafter(1, 0))

# +-------------------------------------------------------------+
# | Computed constants, do not manually modify below this line. |
# +-------------------------------------------------------------+

"""
The base probability of changing our estimate from one family of sequences to
another, part way through adjudicating a sequence.
"""
CHANGE_PROB = 1.0 - SAME_PROB

LOGGED_CHANGE_PROB = math.log10(CHANGE_PROB)

LOGGED_SAME_PROB = math.log10(SAME_PROB)
