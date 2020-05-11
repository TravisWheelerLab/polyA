from typing import Dict, Tuple

# FIXME - GEORGE - this is a collapsed matrix so I updated type to be Dict[Tuple[str, int], int]
# FIXME might break code?
OriginMatrix = Dict[Tuple[str, int], int]

"""
Hash implementation of sparse 2D DP matrix. This is a collapsed matrix. Key is Tuple[str, int] 
with the row(subfam) and column. Holds which cell in previous column the probability in 
the DP matrix came from. Used when doing backtrace through the DP matrix. 

Doesn't have any functions - Populated in FillProbMatrix()
"""
