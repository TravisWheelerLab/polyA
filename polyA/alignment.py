from typing import List, NamedTuple


class Alignment(NamedTuple):
    """
    A container to hold data related to a single alignment.

    >>> a = Alignment("", "ch", 25, 100, 0, 0, 0, 0, 0, [], "", 0, "", 0, 0)
    >>> a.chrom_name
    'ch'
    >>> a.chrom_start
    25
    >>> a.chrom_stop
    100
    >>> a.chrom_length
    76
    """

    subfamily: str
    chrom_name: str
    chrom_start: int
    chrom_stop: int
    score: int
    start: int
    stop: int
    consensus_start: int
    consensus_stop: int
    sequences: List[str]
    strand: str
    flank: int
    sub_matrix_name: str
    gap_init: float
    gap_ext: float

    @property
    def chrom_length(self) -> int:
        return self.chrom_stop - self.chrom_start + 1

    @property
    def sequence(self) -> str:
        return self.sequences[0]

    @sequence.setter
    def sequence(self, value) -> None:
        self.sequences[0] = value

    @property
    def subfamily_sequence(self) -> str:
        return self.sequences[1]

    @subfamily_sequence.setter
    def subfamily_sequence(self, value) -> None:
        self.sequences[1] = value


_skip = Alignment(
    subfamily="skip",
    chrom_name="skip_chrom",
    chrom_start=0,
    chrom_stop=0,
    score=0,
    start=0,
    stop=0,
    consensus_start=0,
    consensus_stop=0,
    sequences=["", ""],
    strand="",
    flank=0,
    sub_matrix_name="",
    gap_init=0,
    gap_ext=0,
)


def get_skip_state():
    """
    Return the skip state singleton.
    """
    return _skip
