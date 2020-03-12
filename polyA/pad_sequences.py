import logging
from typing import Iterable
from .align_score_matrix import Alignment
from .edges import edges


_logger = logging.getLogger(__name__)


def pad_sequences(alignments: Iterable[Alignment]):
    """
    Pad out sequences with "." to allow regions where sequences do not all
    have the same start and stop positions.

    >>> a0 = Alignment(sequences=["a", "b"], start=1, stop=1, subfamily='', score=0, consensus_start=0, consensus_stop=0)
    >>> a1 = Alignment(sequences=["aaa", "bbb"], start=0, stop=2, subfamily='', score=0, consensus_start=0, consensus_stop=0)
    >>> pad_sequences([a0, a1])
    >>> a0.sequence
    '.a.'
    >>> a0.subfamily_sequence
    '.b.'
    >>> a1.sequence
    'aaa'
    >>> a1.subfamily_sequence
    'bbb'
    """
    edgeStart, edgeStop = edges(alignments)

    for alignment in alignments:
        if alignment.sequence == "":
            continue

        leftPad = (alignment.start - edgeStart) * "."
        rightPad = (edgeStop - alignment.stop) * "."

        alignment.sequence = leftPad + alignment.sequence + rightPad
        alignment.subfamily_sequence = (
            leftPad + alignment.subfamily_sequence + rightPad
        )