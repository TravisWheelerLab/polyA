DEFAULT_TRANS_PENALTY = 55
"""
The positive exponent value for the base probability of changing annotations.
A value of 55 means the probability will be 1e-55.
"""

DEFAULT_CHUNK_SIZE = 31
"""
The width of the "window" used to break a sequence up into segments.
Measured in base pairs.
"""

INFINITE_SHARD_GAP = -1
"""
A magic value that indicates no sharding should occur.
"""

DEFAULT_SHARD_GAP = INFINITE_SHARD_GAP
"""
Allowed gap between sequences for them to be included in the same shard.
"""

NAN_STRING = "NaN"
"""
The string used to represent "not a number" for serialization and
deserialization.
"""

PROB_SKIP = 0.4  # about 60% of genome is TE derived
"""
used for prior counts, about 60% of genome is TE derived, so the skip
state will have a prob of 40%
"""

PROB_SKIP_TR = 0.06
"""
About 6% of the genome is expected to be tandem repeats.
"""

SAME_PROB_LOG = 0.0
"""
penalty for staying in the same state in the DP.
mathematically not exactly 0, but in python log(1-10e-45) = 0
so set to 0 and avoid doing the math
"""

SKIP_ALIGN_SCORE = 10.0
"""
alignment score given to skip state, no lambda adjustment needed
"""

START_ID = 1111
"""
The arbitrary temporary ID used to bootstrap the graph process. The
value doesn't matter, but it's here for clarity.
"""

CROSS_MATCH_ADJUSTMENT = 0.999
"""
A constant value given from cross_match that is used to compute 
the complexity adjusted score of an alignment.
"""

# +-------------------------------------------------------------+
# | Computed constants, do not manually modify below this line. |
# +-------------------------------------------------------------+
