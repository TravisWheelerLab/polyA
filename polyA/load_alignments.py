import re
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


def _parse_chrom_meta(line: str) -> Optional[Tuple[str, int, int]]:
    match = re.search(r"(.+):(\d+)-(\d+)", line.strip())
    if match is None:
        return None

    chrom_name = match.groups()[0]
    chrom_start = int(match.groups()[1])
    chrom_stop = int(match.groups()[2])

    return chrom_name, chrom_start, chrom_stop


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
    chrom_name: str = ""
    chrom_start: int = 0
    chrom_stop: int = 0

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

            chrom_meta = _parse_chrom_meta(meta["TR"])
            if chrom_meta is not None:
                chrom_name, chrom_start, chrom_stop = chrom_meta
            else:
                raise ValueError("metadata incomplete, missing TR")

            if meta["TQ"] == "t":
                yield Alignment(
                    subfamily=meta["ID"],
                    chrom_name=chrom_name,
                    chrom_start=chrom_start,
                    chrom_stop=chrom_stop,
                    score=int(meta["SC"]),
                    start=int(meta["ST"]),
                    stop=int(meta["SP"]),
                    consensus_start=int(meta["CST"]),
                    consensus_stop=int(meta["CSP"]),
                    sequences=[seqs[0][::-1], seqs[1][::-1]],
                    strand=meta["SD"],
                    flank=int(meta["FL"]),
                    sub_matrix_name=meta["MX"],
                    gap_init=float(meta["GI"]),
                    gap_ext=float(meta["GE"]),
                )
            else:
                yield Alignment(
                    subfamily=meta["ID"],
                    chrom_name=chrom_name,
                    chrom_start=chrom_start,
                    chrom_stop=chrom_stop,
                    score=int(meta["SC"]),
                    start=int(meta["ST"]),
                    stop=int(meta["SP"]),
                    consensus_start=int(meta["CST"]),
                    consensus_stop=int(meta["CSP"]),
                    sequences=[seqs[0], seqs[1]],
                    strand=meta["SD"],
                    flank=int(meta["FL"]),
                    sub_matrix_name=meta["MX"],
                    gap_init=float(meta["GI"]),
                    gap_ext=float(meta["GE"]),
                )

            meta.clear()
            seqs.clear()
            chrom_name = ""
            chrom_start = 0
            chrom_stop = 0


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
    >>> a0 = Alignment("", "a", 1, 100, 0, 1, 10, 0, 0, [], "", 0, "", 0, 0)
    >>> a1 = Alignment("", "a", 1, 100, 0, 21, 30, 0, 0, [], "", 0, "", 0, 0)
    >>> shards = list(shard_overlapping_alignments([a0, a1], 10))
    >>> len(shards)
    2
    >>> shards[0].start
    1
    >>> shards[0].stop
    15
    >>> len(shards[0].alignments)
    2
    >>> shards[1].start
    16
    >>> shards[1].stop
    100
    >>> len(shards[1].alignments)
    2
    """
    shard_alignments: List[Alignment] = (
        [get_skip_state()] if add_skip_state else []
    )

    shard_start = 1
    shard_stop = None

    is_infinite_gap = shard_gap == INFINITE_SHARD_GAP

    last_alignment = None

    for alignment in alignments:
        last_alignment = alignment

        # If this is the first alignment we need to initialize
        # the stop position
        if shard_stop is None:
            shard_stop = alignment.stop

        if is_infinite_gap or alignment.start <= (shard_stop + shard_gap):
            # This alignment might extend past the previous
            # alignments in the shard, capture that here
            if shard_stop < alignment.stop:
                shard_stop = alignment.stop

            shard_alignments.append(alignment)
        else:
            # The gap between shards might be bigger than the
            # maximum shard gap, so we have to consume half of
            # the actual gap for this shard
            actual_gap = alignment.start - shard_stop
            shard_stop += int(actual_gap / 2)

            yield Shard(
                start=shard_start,
                stop=shard_stop,
                alignments=shard_alignments,
            )

            shard_start = shard_stop + 1
            shard_stop = alignment.stop

            # Note: important to create a new list here or we will
            # mutate the one we just handed back to the caller.
            shard_alignments = (
                [get_skip_state(), alignment] if add_skip_state else [alignment]
            )

    if last_alignment is None:
        raise ValueError("no alignments found")

    yield Shard(
        start=shard_start,
        stop=last_alignment.chrom_stop - last_alignment.chrom_start + 1,
        alignments=shard_alignments,
    )
