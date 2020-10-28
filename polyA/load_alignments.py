from typing import Iterable, List, Optional, TextIO, Tuple
from .alignment import Alignment


def _parse_meta_line(line: str) -> Optional[Tuple[str, str]]:
    """
    Attempt to parse one line of a Stockholm file. If the line
    contains metadata, return a tuple of the key and value.
    Otherwise, return ``None``.

    >>> _parse_meta_line("#=GF ID\tSeqMcSeq")
    ('ID', 'SeqMcSeq')
    >>> _parse_meta_line("seqname\tatcg")
    >>> _parse_meta_line("//")
    >>> _parse_meta_line("# stockholm 1.0")
    """
    if line.strip().startswith("#="):
        parts = line.strip().split()

        if len(parts) != 3:
            return None

        return parts[1], parts[2]

    return None


def _parse_alignment_line(line: str) -> Optional[Tuple[str, str]]:
    """
    Attempt to parse one line of a sequence alignment,
    returning a tuple of the sequence name and the sequence
    itself. If the string isn't formatted correctly, return
    ``None``.

    >>> _parse_alignment_line("seqname\tatcg")
    ('seqname', 'atcg')
    >>> _parse_alignment_line("#=GF ID SeqMcSeq")
    >>> _parse_meta_line("//")
    >>> _parse_meta_line("# stockholm 1.0")
    """
    parts = line.strip().split()

    if len(parts) != 2:
        return None

    return parts[0], parts[1]


def _parse_preamble_line(line: str) -> bool:
    """
    Return ``True`` if this line contains a Stockholm format
    preamble ("# Stockholm 1.0") and ``False`` otherwise.
    """
    return line.strip().upper().startswith("# STOCKHOLM")


def _parse_terminator_line(line: str) -> bool:
    """
    Return ``True`` if this line contains a terminator
    ("//") or ``False`` otherwise.
    """
    return line.strip() == "//"


def load_alignments(file: TextIO) -> Iterable[Alignment]:
    """
    Load a set of alignments in Stockholm format.

    ::

        #=GF ID  MERX#DNA/TcMar-Tigger
        #=GF TR  chr0:0000-0000
        #=GF SC  1153
        #=GF ST  +
        #=GF TQ  -1
        #=GF ST  127
        #=GF SP  601
        #=GF CST 135
        #=GF CSP 628
        #=GF FL  128
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

    meta = {}
    seqs = []

    for line in file:
        if _parse_preamble_line(line):
            continue

        parts = _parse_meta_line(line)
        if parts is not None:
            meta[parts[0]] = parts[1]
            continue

        parts = _parse_alignment_line(line)
        if parts is not None:
            seqs.append(parts[1])
            continue

        if _parse_terminator_line(line):
            # TODO: Validation and good error messages go here
            # Verify that...
            #   we have two sequences
            #   all required metadata was present

            alignments.append(
                Alignment(
                    subfamily=meta["ID"],
                    chrom=meta["TR"],
                    score=int(meta["SC"]),
                    start=int(meta["ST"]),
                    stop=int(meta["SP"]),
                    consensus_start=int(meta["CST"]),
                    consensus_stop=int(meta["CSP"]),
                    sequences=[seqs[0], seqs[1]],
                    strand=meta["SD"],
                    flank=int(meta["FL"]),
                )
            )
            meta.clear()
            seqs.clear()
            continue

    return alignments
