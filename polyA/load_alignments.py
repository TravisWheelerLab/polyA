from typing import Iterable, List, NamedTuple, Optional, TextIO, Tuple

from .alignment import Alignment, get_skip_state
from .constants import INFINITE_SHARD_GAP


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

    >>> _parse_terminator_line(" // ")
    True
    >>> _parse_terminator_line(" / / ")
    False
    """
    return line.strip() == "//"


def load_alignments(
    file: TextIO, add_skip_state: bool = False
) -> Iterable[Alignment]:
    """
    Load a set of alignments in Stockholm format.
    See below for an example.

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
        #=GF MX  20p41g.matrix
        #=GF GI -25
        #=GF GE -5
    """
    if add_skip_state:
        yield get_skip_state()

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

            if meta["TQ"] == "t":
                yield Alignment(
                    subfamily=meta["ID"],
                    chrom=meta["TR"],
                    score=int(meta["SC"]),
                    start=int(meta["ST"]),
                    stop=int(meta["SP"]),
                    consensus_start=int(meta["CST"]),
                    consensus_stop=int(meta["CSP"]),
                    sequences=[seqs[0][::-1], seqs[1][::-1]],
                    strand=meta["SD"],
                    flank=int(meta["FL"]),
                    sub_matrix_name=meta["MX"],
                    gap_init=int(meta["GI"]),
                    gap_ext=int(meta["GE"]),
                )
            else:
                yield Alignment(
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
                    sub_matrix_name=meta["MX"],
                    gap_init=int(meta["GI"]),
                    gap_ext=int(meta["GE"]),
                )
            meta.clear()
            seqs.clear()


class Shard(NamedTuple):
    start: int
    stop: int
    alignments: List[Alignment]


def shard_overlapping_alignments(
    alignments: Iterable[Alignment],
    shard_gap: int,
    add_skip_state: bool = True,
) -> Iterable[Shard]:
    """
    Shard the given alignments into overlapping groups, separated by no more
    than `shard_gap` nucleotides. This allows for more efficient processing
    for alignments over large regions of sequence (such as an entire genome)
    where many regions will be empty.
    Precondition: alignments are sorted by their start position and all
    alignments have start position <= stop position.
    >>> skip = get_skip_state()
    >>> a0 = Alignment("", "", 0, 0, 10, 0, 0, [], "", 0, "", 0, 0)
    >>> a1 = Alignment("", "", 0, 2, 12, 0, 0, [], "", 0, "", 0, 0)
    >>> a2 = Alignment("", "", 0, 12 + 50, 70, 0, 0, [], "", 0, "", 0, 0)
    >>> a3 = Alignment("", "", 0, 70 + 51, 140, 0, 0, [], "", 0, "", 0, 0)
    >>> a4 = Alignment("", "", 0, 70 + 51, 140, 0, 0, [], "", 0, "", 0, 0)
    >>> a5 = Alignment("", "", 0, 70 + 52, 130, 0, 0, [], "", 0, "", 0, 0)
    >>> chunks = list(shard_overlapping_alignments([a0, a1, a2, a3, a4, a5], 50))
    >>> len(chunks)
    2
    >>> chunks[0].start
    0
    >>> chunks[0].stop
    95
    >>> chunks[0].alignments == [skip, a0, a1, a2]
    True
    >>> chunks[1].start
    96
    >>> chunks[1].stop
    140
    >>> chunks[1].alignments == [skip, a3, a4, a5]
    True
    """
    shard_alignments: List[Alignment] = (
        [get_skip_state()] if add_skip_state else []
    )
    shard_start = 0
    shard_stop = None

    for alignment in alignments:
        is_start = shard_stop is None
        is_infinite_gap = shard_gap == INFINITE_SHARD_GAP

        if is_start:
            shard_stop = alignment.stop

        if (
            is_start
            or is_infinite_gap
            or alignment.start <= (shard_stop + shard_gap)
        ):
            if shard_stop < alignment.stop:
                shard_stop = alignment.stop

            shard_alignments.append(alignment)
        else:
            shard_stop += int(shard_gap / 2)
            yield Shard(
                start=shard_start,
                stop=shard_stop,
                alignments=shard_alignments,
            )

            shard_start = shard_stop + 1
            shard_stop = shard_start

            # Note: important to create a new list here or we will
            # mutate the one we just handed back to the caller.
            shard_alignments = (
                [get_skip_state(), alignment] if add_skip_state else [alignment]
            )

    yield Shard(
        start=shard_start,
        stop=shard_alignments[-1].stop,
        alignments=shard_alignments,
    )


def shard_overlapping_alignments_update(
    alignments: Iterable[Alignment],
    shard_gap: int,
    add_skip_state: bool = True,
) -> Iterable[Shard]:
    """
    Shard alignments with bug fix

    """
    shard_alignments: List[Alignment] = (
        [get_skip_state()] if add_skip_state else []
    )
    shard_start = 0
    shard_stop = None  # None did not work here because of line 215

    for alignment in alignments:
        is_start = shard_stop is None
        is_infinite_gap = shard_gap == INFINITE_SHARD_GAP
        if (
            is_start
            or is_infinite_gap
            or alignment.start <= (shard_stop + shard_gap)
        ):
            if is_start or shard_stop < alignment.stop:
                # print("change")
                shard_stop = alignment.stop

            shard_alignments.append(alignment)
        else:
            shard_stop += int(shard_gap / 2)
            yield Shard(
                start=shard_start,
                stop=shard_stop,
                alignments=shard_alignments,
            )

            shard_start = shard_stop + 1
            shard_stop = None

            # Note: important to create a new list here or we will
            # mutate the one we just handed back to the caller.
            shard_alignments = (
                [get_skip_state(), alignment] if add_skip_state else [alignment]
            )
    yield Shard(
        start=shard_start,
        stop=shard_alignments[-1].stop,
        alignments=shard_alignments,
    )
