from polyA.fill_align_matrix import fill_align_matrix
from polyA.edges import edges
from polyA.constants import SKIP_ALIGN_SCORE
from polyA.lambda_provider import ConstantLambdaProvider
from polyA.substitution_matrix import load_substitution_matrices, SubMatrix


def test_fill_align_matrix():
    chrom_seqs = [
        "",
        "TCAGACTGTTCATGAGTGCTCACCTGGTAGAGG-----AAAGACCCTG-AAC----GTCCAGAGCTTTCCCCTGACTGGCCACGTGCGGACTCTTTGTGGCT---CTGCTGA................",
        "CCATCTTTAAGAATGGAAGATGGACATGAG-AAACTTTTCT--CTA-----CCATCTACATAAAACAAATGA--ATAGAGTTTTCCTGACTTAA---ATAAAATTGCCCATTTTTTGGCCTTGTGCC-CACTTA------CAGCTGGACTTAGTAGTCTCACA--GACCTCAGTGTGCTCTAGGTAAGATCATCCCTCTTGCTGATCTTGCCATTGCTAGCAACAAGAATTTAGAAAGGCTGACT................",
        "TTAAATGACTTTTAATGAAGGGGTTT-TCGATTGAATTATTTACTAATTAAACGA--TTCTGTGTTAAGAGACTTTCGGTCAATCGGTATACTGTTCTATGAGTTACT-TTCAAG-----TAA--GTTAAAAAA----TAAACAA--AAAATTTGTAAAGTACTTCAAATATTATAAAG--TAACATATCCTACCTAAATTATCGAGAAGTTATAAATTTCTCTTTCCAAAA-AGACTGACCCAAAGTAAAAC................",
    ]
    subfam_seqs = [
        "",
        "TCAGACTGTTCA-----ACTCACCTGGCAGCCACTTCCAGAGCCCCTGGAACTCTAGCCCAAGGCTCTC---TGACTGACCCCTTCTGAGATCTTCTTGGCTTAGCAGCTGA................",
        "CCAGCCTTAGGACTGCCAGAT----ATAACTAAGCCTTTCTTTCTATATGTCCATGAACG-AAAAGGAATGGCTATAGGGGT--CCTGACTCAAGTCATAGGAT----CATGGCATCGCCCGGTGGCATACTTATCCAAGCAACTGGACTCCGTGGCNCTAGGATGGCCTC--CTTGCCTTAGG---G--CA-----CTAGCTGNCACTGCC-CTACTGGCA-CAAGAAGCTAACAAA-CTGACT................",
        "TTAAGCGGCTTTCGGTTAAGCGGCTTACCGGTTAAGCGGCTT-TTAGTTAAGCGGCTTTCTG-GTTAAGCGGCTTTCGGTTAAGCGGCTTTTGGTTAAGCGGCTTACTGGTTAAGCGGTTTAATGGTTAAATGGTTTTTGAACAAAGAAAATTGATAAA-TATTTCAAATATCGTTAAGCATAACCTACCCTACCAAAATTGTCGGAAACTCCTAAA-AACGCTT--CAAAATCAACTGGCCAAAANTAAAAC................",
    ]
    alignment_start_positions = [6409, 6410, 6450, 6463]
    alignment_stop_positions = [6696, 6508, 6672, 6695]
    gap_inits = [0, -25.0, -25.0, -25.0]
    gap_exts = [0, -5.0, -5.0, -5.0]
    start_all, stop_all = edges(
        alignment_start_positions, alignment_stop_positions
    )
    column_count = stop_all - start_all + 1 + 2

    _lambda_provider = ConstantLambdaProvider(0.1227)
    with open(
        f"fixtures/ultra_test_files/ex13.fa.cm.matrix", "r"
    ) as _sub_matrices_file:
        sub_matrices = load_substitution_matrices(
            _sub_matrices_file, _lambda_provider, False
        )
    lambda_values = [0.0]
    alignment_substitution_matrices = [SubMatrix("skip", 0.0)]

    for seq in chrom_seqs:
        lambda_values.append(sub_matrices["matrix1"].lamb)
        alignment_substitution_matrices.append(sub_matrices["matrix1"])

    align_matrix = fill_align_matrix(
        lambda_values,
        column_count,
        start_all,
        31,
        gap_inits,
        gap_exts,
        SKIP_ALIGN_SCORE,
        subfam_seqs,
        chrom_seqs,
        alignment_start_positions,
        alignment_stop_positions,
        [sm.scores for sm in alignment_substitution_matrices],
        [sm.background_freqs for sm in alignment_substitution_matrices],
    )
    # check values of align matrix
    # assert certain values exists or don't exist
    assert (1, 0) not in align_matrix
    assert (1, 100) not in align_matrix
    assert align_matrix[1, 1] == 15.690262500000001
    assert align_matrix[1, 99] == 12.59975625

    assert (2, 40) not in align_matrix
    assert (2, 264) not in align_matrix
    assert align_matrix[2, 41] == 13.075218750000001
    assert align_matrix[2, 263] == 4.99235625

    assert (3, 53) not in align_matrix
    assert (3, 287) not in align_matrix
    assert align_matrix[3, 54] == 17.11665
    assert align_matrix[3, 286] == 22.822200000000002

    assert align_matrix[0, 0] == 10.0
    assert align_matrix[0, 287] == 10.0


def test_fill_align_matrix_complexity_adjusted():
    subfam_seqs = [
        "",
        "TCAGACTGTTCA-----ACTCACCTGGCAGCCACTTCCAGAGCCCCTGGAACTCTAGCCCAAGGCTCTC---TGACTGACCCCTTCTGAGATCTTCTTGGCTTAGCAGCTGA................",
        "CCAGCCTTAGGACTGCCAGAT----ATAACTAAGCCTTTCTTTCTATATGTCCATGAACG-AAAAGGAATGGCTATAGGGGT--CCTGACTCAAGTCATAGGAT----CATGGCATCGCCCGGTGGCATACTTATCCAAGCAACTGGACTCCGTGGCNCTAGGATGGCCTC--CTTGCCTTAGG---G--CA-----CTAGCTGNCACTGCC-CTACTGGCA-CAAGAAGCTAACAAA-CTGACT................",
        "TTAAGCGGCTTTCGGTTAAGCGGCTTACCGGTTAAGCGGCTT-TTAGTTAAGCGGCTTTCTG-GTTAAGCGGCTTTCGGTTAAGCGGCTTTTGGTTAAGCGGCTTACTGGTTAAGCGGTTTAATGGTTAAATGGTTTTTGAACAAAGAAAATTGATAAA-TATTTCAAATATCGTTAAGCATAACCTACCCTACCAAAATTGTCGGAAACTCCTAAA-AACGCTT--CAAAATCAACTGGCCAAAANTAAAAC................",
    ]
    chrom_seqs = [
        "",
        "TCAGACTGTTCATGAGTGCTCACCTGGTAGAGG-----AAAGACCCTG-AAC----GTCCAGAGCTTTCCCCTGACTGGCCACGTGCGGACTCTTTGTGGCT---CTGCTGA................",
        "CCATCTTTAAGAATGGAAGATGGACATGAG-AAACTTTTCT--CTA-----CCATCTACATAAAACAAATGA--ATAGAGTTTTCCTGACTTAA---ATAAAATTGCCCATTTTTTGGCCTTGTGCC-CACTTA------CAGCTGGACTTAGTAGTCTCACA--GACCTCAGTGTGCTCTAGGTAAGATCATCCCTCTTGCTGATCTTGCCATTGCTAGCAACAAGAATTTAGAAAGGCTGACT................",
        "TTAAATGACTTTTAATGAAGGGGTTT-TCGATTGAATTATTTACTAATTAAACGA--TTCTGTGTTAAGAGACTTTCGGTCAATCGGTATACTGTTCTATGAGTTACT-TTCAAG-----TAA--GTTAAAAAA----TAAACAA--AAAATTTGTAAAGTACTTCAAATATTATAAAG--TAACATATCCTACCTAAATTATCGAGAAGTTATAAATTTCTCTTTCCAAAA-AGACTGACCCAAAGTAAAAC................",
    ]
    alignment_start_positions = [6409, 6410, 6450, 6463]
    alignment_stop_positions = [6696, 6508, 6672, 6695]
    gap_inits = [0, -25.0, -25.0, -25.0]
    gap_exts = [0, -5.0, -5.0, -5.0]
    start_all, stop_all = edges(
        alignment_start_positions, alignment_stop_positions
    )
    column_count = stop_all - start_all + 1 + 2

    _lambda_provider = ConstantLambdaProvider(0.1227)
    with open(
        f"fixtures/ultra_test_files/ex13.fa.cm.matrix", "r"
    ) as _sub_matrices_file:
        sub_matrices = load_substitution_matrices(
            _sub_matrices_file, _lambda_provider, True
        )
    lambda_values = [0.0]
    alignment_substitution_matrices = [SubMatrix("skip", 0.0)]

    for seq in chrom_seqs:
        lambda_values.append(sub_matrices["matrix1"].lamb)
        alignment_substitution_matrices.append(sub_matrices["matrix1"])

    align_matrix = fill_align_matrix(
        lambda_values,
        column_count,
        start_all,
        31,
        gap_inits,
        gap_exts,
        SKIP_ALIGN_SCORE,
        subfam_seqs,
        chrom_seqs,
        alignment_start_positions,
        alignment_stop_positions,
        [sm.scores for sm in alignment_substitution_matrices],
        [sm.background_freqs for sm in alignment_substitution_matrices],
    )

    assert (1, 0) not in align_matrix
    assert (1, 100) not in align_matrix
    assert align_matrix[1, 1] == 14.762038466561917
    assert align_matrix[1, 99] == 11.333849987455658

    assert (2, 40) not in align_matrix
    assert (2, 264) not in align_matrix
    assert align_matrix[2, 41] == 12.98841727761375
    assert align_matrix[2, 263] == 4.912091167445111

    assert (3, 53) not in align_matrix
    assert (3, 287) not in align_matrix
    assert align_matrix[3, 54] == 15.686524697939982
    assert align_matrix[3, 286] == 21.288106662982315

    assert align_matrix[0, 0] == 10.0
    assert align_matrix[0, 287] == 10.0
