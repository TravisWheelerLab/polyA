from __future__ import annotations
from itertools import groupby
from typing import Callable, Iterable, List, NamedTuple, TextIO, Union


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
    sequence: str
    subfamily_sequence: str

    def update(
        self, sequence: Union[str, None], subfamily_sequence: Union[str, None]
    ) -> Alignment:
        return Alignment(
            subfamily=self.subfamily,
            score=self.score,
            start=self.start,
            stop=self.stop,
            consensus_start=self.consensus_start,
            consensus_stop=self.consensus_stop,
            sequence=self.sequence if sequence is None else sequence,
            subfamily_sequence=self.subfamily_sequence
            if subfamily_sequence is None
            else subfamily_sequence,
        )


def _line_grouper(prefix: str) -> Callable[[str], str]:
    """
    Helper function to return a function appropriate for use with
    `itertools.groupby` that checks for the key prefix and includes all keys
    themselves under a key associated with the empty string.
    """
    lastKey = ""

    def _group_lines(line: str) -> str:
        nonlocal lastKey
        if line.startswith(prefix):
            lastKey = line
            return ""
        return lastKey

    return _group_lines


def load_alignments(file: TextIO) -> Iterable[Alignment]:
    """
    Load a set of alignments from a file formatted the way cross_match
    formats its output (see http://www.phrap.org/phredphrapconsed.html).
    For example:

    Align: AluYj4	1-AluJr_311-AluYb8_629-AluSz6	1049	+	1	309	1	311
    >1-AluJr_311-AluYb8_629-AluSz6
    GGCCTTGCGAGGTGGGTCACGNCANNTGTAATCCCACTAATTTGGCCGGCCGAGGGTGGC
    GGATC----GTTCANNNAGATTTTGAGGCCAGCCTGGGGGACCTANNN--GCGAGAGGCC
    ...
    >AluYj4
    GGCCGGGCGCGGTGGCTCGCGCC---TGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGC
    GGATCACGAGGTCAGG-AGATC--GAGACCATCCTGG-----CTAACACGGTGAAACCCC
    ...
    """
    for meta, lines in groupby(file, _line_grouper("Align:")):
        if meta == "":
            continue
        sequences: List[str] = []
        for element, lines in groupby(lines, _line_grouper(">")):
            if element == "":
                continue
            sequence = "".join(map(lambda l: l.strip(), lines))
            sequences.append(sequence)

        metaItems = meta.split()

        subfamily = metaItems[1]
        score = int(metaItems[3])
        strand = metaItems[4]

        start = int(metaItems[5])
        stop = int(metaItems[6])
        if start > stop:
            start, stop = stop, start

        consensusStart = int(metaItems[7])
        consensusStop = int(metaItems[8])
        if consensusStart > consensusStop:
            consensusStart, consensusStop = consensusStop, consensusStart

        sequence = sequences[0]
        subfamilySequence = sequences[1]
        if strand == "-":
            sequence = sequence[-1::-1]
            subfamilySequence = subfamilySequence[-1::-1]

        yield Alignment(
            subfamily=subfamily,
            score=score,
            start=start,
            stop=stop,
            consensus_start=consensusStart,
            consensus_stop=consensusStop,
            sequence=sequence,
            subfamily_sequence=subfamilySequence,
        )
