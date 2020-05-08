from typing import Iterable, List, NamedTuple


class Alignment(NamedTuple):
    """
    A container to hold data related to a single alignment.
    """

    subfamily: str
    score: int
    start: int
    stop: int
    consensus_start: int
    consensus_stop: int
    sequences: List[str]
    strand: str

    def __repr__(self):
        return f"Alignment(subfamily={self.subfamily})"

    def __str__(self):
        return f"Alignment(subfamily={self.subfamily})"

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
