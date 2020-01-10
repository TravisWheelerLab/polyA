from typing import Any, Callable, Dict, List, Tuple

"""
Typedef to represent a sparse confidence matrix implemented as a
dictionary (for now).
"""
ConfScoreMatrix = Dict[Tuple[int, int], float]
