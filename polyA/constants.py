import numpy as np

DEFAULT_CHUNK_SIZE = 30
"""
The width of the "window" used to break a sequence up into pieces.
Measured in base pairs.
"""

DEFAULT_GAP_EXTEND = -5
"""
TODO: Explanation (Kaitlin)
"""

DEFAULT_GAP_START = -25
"""
TODO: Explanation (Kaitlin)
"""

DEFAULT_LAMBDA = 0.1227
"""
TODO: Explanation (Kaitlin)
TODO: We don't want a default here, but then we require Easel
"""

NAN_STRING = "NaN"
"""
The string used to represent "not a number" for serialization and
deserialization.
"""

DEFAULT_CHANGE_PROB = float(np.nextafter(0, 1))
"""
The base probability that we change our estimate from one
base pair to the next. This will generally be a very small value.
"""

# +-------------------------------------------------------------+
# | Computed constants, do not manually modify below this line. |
# +-------------------------------------------------------------+
