import math
import sys
import numpy as np

"""
The width of the "window" used to break a sequence up into pieces.
Measured in base pairs.
"""
DEFAULT_CHUNK_SIZE = 30

"""
TODO: Explanation (Kaitlin)
"""
DEFAULT_GAP_EXT = -5

"""
TODO: Explanation (Kaitlin)
"""
DEFAULT_GAP_INIT = -25

"""
The string used to represent "not a number" for serialization and
deserialization.
"""
NAN_STRING = "NaN"

"""
The base probability that we do not change our estimate from one
base pair to the next. This will generally be a very small value.
"""
DEFAULT_SAME_PROB = float(np.nextafter(1, 0))

# +-------------------------------------------------------------+
# | Computed constants, do not manually modify below this line. |
# +-------------------------------------------------------------+

"""
The base probability of changing our estimate from one family of
sequences to another, part way through adjudicating a sequence.
"""
DEFAULT_CHANGE_PROB = 1.0 - DEFAULT_SAME_PROB
