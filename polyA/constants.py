DEFAULT_CHUNK_SIZE = 31
"""
The width of the "window" used to break a sequence up into pieces.
Measured in base pairs.
"""

DEFAULT_GAP_EXTEND = -5
"""
Penalty for extending a gap in the alignment
"""

DEFAULT_GAP_START = -25
"""
Penalty for opening a gap in the alignment
"""

DEFAULT_LAMBDA = 0.1227
"""
scaling factor that corresponds to the substitution matrix used, needed for confidence 
calculations

TODO: We don't want a default here, but then we require Easel
TODO: Instead of a default here we want it to be input by the user at the command 
line - lambda corresponds to the substitution matrix used when scoring alignments - 
then if the user doesn't input anything we use easel to calculate it
"""

NAN_STRING = "NaN"
"""
The string used to represent "not a number" for serialization and
deserialization.
"""

DEFAULT_CHANGE_PROB = 10**-45
"""
The probability that corresponds to the penalty of changing rows (jumping subfams)
in the DP matrix. This will generally be a very small value.
"""

# +-------------------------------------------------------------+
# | Computed constants, do not manually modify below this line. |
# +-------------------------------------------------------------+
