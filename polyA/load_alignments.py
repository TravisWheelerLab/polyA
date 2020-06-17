from itertools import groupby
from typing import Callable, Iterable, List, TextIO
from .alignment import Alignment
from .pad_sequences import pad_sequences


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
    For example::

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
    alignments: List[Alignment] = [
        Alignment(
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
        )
    ]
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
        chrom = metaItems[2]
        score = int(metaItems[3])

        # TODO: Figure out whether we can drop this bit
        strand = metaItems[4]

        start = int(metaItems[5])
        stop = int(metaItems[6])

        consensusStart = int(metaItems[7])
        consensusStop = int(metaItems[8])
        flank = int(metaItems[9])

        sequence = sequences[0]
        subfamilySequence = sequences[1]

        alignments.append(
            Alignment(
                subfamily=subfamily,
                chrom=chrom,
                score=score,
                start=start,
                stop=stop,
                consensus_start=consensusStart,
                consensus_stop=consensusStop,
                sequences=[sequence, subfamilySequence],
                strand=strand,
                flank=flank,
            )
        )

    # pad_sequences(alignments)

    return alignments
