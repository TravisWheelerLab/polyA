CHANGE_PROB = 10 ** -45
"""
The base probability of changing annotations.
"""

DEFAULT_CHUNK_SIZE = 31
"""
The width of the "window" used to break a sequence up into segments.
Measured in base pairs.
"""

DEFAULT_GAP_INIT = -25
"""
Penalty given to start a gap in alignment
"""

DEFAULT_GAP_EXT = -5
"""
Penalty given to extend a gap in alignment
"""

NAN_STRING = "NaN"
"""
The string used to represent "not a number" for serialization and
deserialization.
"""

PROB_SKIP = 0.4  # about 60% of genome is TE derived
"""
TODO(Kaitlin) Explain this better and rename if necessary
"""

PROB_SKIP_TR = 0.06
"""
About 6% of the genome is expected to be tandem repeats.
"""

SAME_PROB_LOG = 0.0
"""
TODO(Kaitlin): Explain this and why we set it to zero
"""

SKIP_ALIGN_SCORE = 5.0
"""
TODO(Kaitlin): Explain this and why we use this value
"""

START_ID = 1111
"""
The arbitrary temporary ID used to bootstrap the graph process. The
value doesn't matter, but it's here for clarity.
"""

# +-------------------------------------------------------------+
# | Computed constants, do not manually modify below this line. |
# +-------------------------------------------------------------+
