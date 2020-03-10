from typing import Iterable, Tuple
from .alignment import Alignment


def edges(alignments: Iterable[Alignment]) -> Tuple[int, int]:
    """
    Find and return the min and max stop positions for the entire region
    included in the alignments.

    >>> a0 = Alignment("sub", 0, 3, 4, 3, 4, ["foobar", "foobaz"])
    >>> a1 = Alignment("sub", 0, 2, 3, 2, 3, ["foobar", "foobaz"])
    >>> a2 = Alignment("sub", 0, 1, 2, 1, 2, ["", ""])
    >>> n, x = edges((a0, a1, a2))
    >>> n
    2
    >>> x
    4
    """
    return (
        min([a.start for a in alignments if a.sequence != ""]),
        max([a.stop for a in alignments if a.sequence != ""]),
    )
