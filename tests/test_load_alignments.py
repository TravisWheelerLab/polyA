from polyA.load_alignments import load_alignments
from pytest import mark


def test_load_alignments_1():
    with open(f"fixtures/alignment1.align", "r") as file:
        skip, first, second = load_alignments(file)

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


def test_load_alignments_2():
    with open(f"fixtures/alignment2.align", "r") as file:
        skip, first, second = load_alignments(file)

    assert first.sequence.startswith("AACAAAAANN")
    assert second.sequence.startswith("AACAAAAANN")

    assert first.subfamily_sequence.startswith("AAAAAAAAA")
    assert second.subfamily_sequence.startswith("AAAAAAAAA")
