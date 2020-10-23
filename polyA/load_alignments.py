from itertools import groupby
from typing import Callable, Iterable, List, TextIO
from .alignment import Alignment


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
        ),
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

        strand = metaItems[4]

        # if on reverse strand, which seq is reversed
        reverse = metaItems[5]

        start = int(metaItems[6])
        stop = int(metaItems[7])

        consensusStart = int(metaItems[8])
        consensusStop = int(metaItems[9])
        flank = int(metaItems[10])

        if reverse == "t":
            sequence = sequences[0][::-1]
            subfamilySequence = sequences[1][::-1]
        else:
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
