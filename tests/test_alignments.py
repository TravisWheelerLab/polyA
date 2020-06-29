from polyA import load_alignments
from pytest import mark

# TODO: Add test files that use different strand values
# TODO: Add test files where start > stop
# TODO: Add test files where consensus_start > consensus_stop


@mark.prob_matrix
def test_alignment_1():
    with open(f"fixtures/alignment1.align", "r") as file:
        first, second = list(load_alignments(file))

    assert first.subfamily == "AluYj4"
    assert second.subfamily == "AluYi6_4d"

    assert first.sequence.startswith("GGCCTTGCGA")
    assert second.sequence.startswith("GGCCTTGCGA")

    assert first.subfamily_sequence.startswith("GGCCGGGCG")
    assert second.subfamily_sequence.startswith("GGCCGGGCG")

    assert first.score == 1049
    assert second.score == 993

    assert first.start == 1
    assert second.start == 1

    assert first.stop == 309
    assert second.stop == 309

    assert first.consensus_start == 1
    assert second.consensus_start == 1

    assert first.consensus_stop == 311
    assert second.consensus_stop == 311

    assert first.flank == 0
    assert second.flank == 0
