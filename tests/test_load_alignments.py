from polyA.alignment import Alignment, get_skip_state
from polyA.load_alignments import load_alignments, shard_overlapping_alignments


def test_load_alignments_1():
    with open(f"fixtures/alignment1.sto", "r") as file:
        skip, first, second = load_alignments(file, add_skip_state=True)

    assert first.subfamily == "L2#LINE/L2"
    assert second.subfamily == "Charlie12#DNA/hAT-Charlie"

    assert first.sequence.startswith("ATTTCCAG---")
    assert second.sequence.startswith("AAGCTTGTCCAA")

    assert first.subfamily_sequence.startswith("ATCTCCAGCCCA")
    assert second.subfamily_sequence.startswith("AAGCTTGTCCA")

    assert first.score == 245
    assert second.score == 671

    assert first.start == 323
    assert second.start == 983

    assert first.stop == 499
    assert second.stop == 1095

    assert first.consensus_start == 2349
    assert second.consensus_start == 2

    assert first.consensus_stop == 2534
    assert second.consensus_stop == 118

    assert first.flank == 548
    assert second.flank == 2755

    assert first.sub_matrix_name == "matrix1"
    assert second.sub_matrix_name == "matrix2"

    assert first.gap_init == 1
    assert second.gap_init == 2

    assert first.gap_ext == 11
    assert second.gap_ext == 22


# test to make sure seqs get flipped with reverse is on target (TQ == 't')
def test_load_alignments_2():
    with open(f"fixtures/alignment2.sto", "r") as file:
        skip, first = load_alignments(file, add_skip_state=True)

    assert first.sequence.startswith("AACAAGAA")
    assert first.subfamily_sequence.startswith("AAAAAAAAA")


def test_load_alignments_no_skip_state():
    with open(f"fixtures/alignment2.sto", "r") as file:
        alignments = list(load_alignments(file, add_skip_state=False))

    assert alignments[0] is not get_skip_state()


def test_shard_alignments_simple():
    # The first shard looks like this:
    #   1   5    10   15   20
    # 0 |--------|
    # 1     |---------|
    # 2  |------|
    s00 = Alignment("", "a", 1, 1000, 0, 1, 10, 0, 0, [], "", 0, "", 0, 0)
    s01 = Alignment("", "a", 1, 1000, 0, 5, 15, 0, 0, [], "", 0, "", 0, 0)
    s02 = Alignment("", "a", 1, 1000, 0, 2, 8, 0, 0, [], "", 0, "", 0, 0)

    # The second shard looks the same, but it is shifted
    # to the right by 115 positions (100 positions beyond
    # right-most position included in the first shard).
    offset = 100 + 15
    s10 = Alignment(
        "",
        "a",
        1,
        1000,
        0,
        s00.start + offset,
        s00.stop + offset,
        0,
        0,
        [],
        "",
        0,
        "",
        0,
        0,
    )
    s11 = Alignment(
        "",
        "a",
        1,
        1000,
        0,
        s01.start + offset,
        s01.stop + offset,
        0,
        0,
        [],
        "",
        0,
        "",
        0,
        0,
    )
    s12 = Alignment(
        "",
        "a",
        1,
        1000,
        0,
        s02.start + offset,
        s02.stop + offset,
        0,
        0,
        [],
        "",
        0,
        "",
        0,
        0,
    )

    skip = get_skip_state()

    shards = list(
        shard_overlapping_alignments([s00, s01, s02, s10, s11, s12], 50)
    )

    assert len(shards) == 2

    assert shards[0].start == 1
    assert shards[0].stop == 65
    assert shards[0].alignments == [skip, s00, s01, s02]

    assert shards[1].start == 66
    assert shards[1].stop == 1000
    assert shards[1].alignments == [skip, s10, s11, s12]


def test_shard_alignments_fixture():
    with open(f"fixtures/ultra_test_files/ex27.fa.cm.sto") as file:
        alignments = load_alignments(file, add_skip_state=False)
        shards = list(shard_overlapping_alignments(alignments, 50, False))

    shard0 = shards[0]
    assert len(shard0.alignments) == 12

    shard1 = shards[1]
    assert len(shard1.alignments) == 27

    shard_last = shards[-1]
    assert shard_last.stop == 146336591 - 146331457 + 1
