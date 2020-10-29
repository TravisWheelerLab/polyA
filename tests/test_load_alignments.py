from polyA.load_alignments import load_alignments
from pytest import mark


def test_load_alignments_1():
    with open(f"fixtures/alignment1.sto", "r") as file:
        skip, first, second = load_alignments(file)

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


#test to make sure seqs get flipped with reverse is on target (TQ == 't')
def test_load_alignments_2():
    with open(f"fixtures/alignment2.sto", "r") as file:
        skip, first = load_alignments(file)

    assert first.sequence.startswith("AACAAGAA")
    assert first.subfamily_sequence.startswith("AAAAAAAAA")

