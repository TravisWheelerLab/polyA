import re
from typing import Dict, List, NamedTuple


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
    sub_matrix_name: str
    gap_init: int
    gap_ext: int
    kimura_divergence: float

    chrom_meta: Dict[str, str] = {}

    # TODO: George - Remove these mutators once we clean up the padding code

    def _read_chrom_meta(self):
        if len(self.chrom_meta) > 0:
            return

        match = re.search(r"(.+):(\d+)-(\d+)", self.chrom)
        self.chrom_meta["name"] = match.groups()[0]
        self.chrom_meta["start"] = match.groups()[1]
        self.chrom_meta["end"] = match.groups()[2]

    @property
    def chrom_name(self) -> str:
        self._read_chrom_meta()
        return self.chrom_meta["name"]

    @property
    def chrom_start(self) -> int:
        self._read_chrom_meta()
        return int(self.chrom_meta["start"])

    @property
    def chrom_end(self) -> int:
        self._read_chrom_meta()
        return int(self.chrom_meta["end"])

    @property
    def chrom_length(self) -> int:
        return self.chrom_end - self.chrom_start

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
    chrom="",
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
    kimura_divergence=0,
)


def get_skip_state():
    """
    Return the skip state singleton.
    """
    return _skip
