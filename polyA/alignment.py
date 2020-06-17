from typing import Callable, Iterable, List, NamedTuple, TextIO


class Alignment(NamedTuple):
    """
    A container to hold data related to a single alignment.
    """

    subfamily: str
    chrom: str
    score: int
    start: int
    stop: int
    consensus_start: int
    consensus_stop: int
    sequences: List[str]
    strand: str
    flank: int

    # TODO: Remove these mutators once we clean up the padding code

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
