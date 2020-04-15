from itertools import groupby
from typing import Callable, Iterable, List, TextIO
from .alignment import Alignment
from .pad_sequences import pad_sequences


def _line_grouper(prefix: str) -> Callable[[str], str]:
    """
    Helper function to return a function appropriate for use with
    `itertools.groupby` that checks for the key prefix and includes all keys
    themselves under a key associated with the empty string.

    >>> grouper = _line_grouper(">")
    >>> grouper(">foo")
    ''
    >>> grouper("bar")
    '>foo'
    >>> grouper("baz")
    '>foo'
    >>> grouper(">buz")
    ''
    >>> grouper("biz")
    '>buz'
    """
    last_key = ""

    def _group_lines(line: str) -> str:
        nonlocal last_key
        if line.startswith(prefix):
            last_key = line
            return ""
        return last_key

    return _group_lines


# TODO: Move implementation to a private function and accept a format param
def load_alignments(file: TextIO) -> Iterable[Alignment]:
    """
    Load a set of pair-wise alignments from a file formatted the way
    cross_match formats its output
    (see http://www.phrap.org/phredphrapconsed.html). For example:

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
            score=0,
            start=0,
            stop=0,
            consensus_start=0,
            consensus_stop=0,
            sequences=["", ""],
            strand="",
        )
    ]
    for meta, meta_lines in groupby(file, _line_grouper("Align:")):
        if meta == "":
            continue
        sequences: List[str] = []
        for element, element_lines in groupby(meta_lines, _line_grouper(">")):
            # Since the alignments are pairwise we'll always go through
            # this loop twice, first for the sequence and second for the
            # subfamily sequence.
            if element == "":
                continue
            sequence = "".join(map(lambda l: l.strip(), element_lines))
            sequences.append(sequence)

        meta_items = meta.split()

        subfamily = meta_items[1]
        score = int(meta_items[3])

        # TODO: Figure out whether we can drop this bit
        strand = meta_items[4]

        start = int(meta_items[5])
        stop = int(meta_items[6])

        consensus_start = int(meta_items[7])
        consensus_stop = int(meta_items[8])

        sequence = sequences[0]
        subfamily_sequence = sequences[1]

        alignments.append(
            Alignment(
                subfamily=subfamily,
                score=score,
                start=start,
                stop=stop,
                consensus_start=consensus_start,
                consensus_stop=consensus_stop,
                sequences=[sequence, subfamily_sequence],
                strand=strand,
            )
        )

    # pad_sequences(alignments)

    return alignments
