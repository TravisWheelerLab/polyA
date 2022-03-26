from typing import Dict

from .alignment import Alignment


def remove_cg_scores(a: Alignment, sub_matrix_scoring: Dict[str, int]) -> int:
    # remove CG score from a.score and return
    remove_score = 0
    consensus = a.sequences[1]
    target = a.sequences[0]
    print(consensus)
    i = 0
    while i < len(consensus) - 1:
        offset = 2
        # check if it contains 2 chars
        char_window = consensus[i : i + offset]
        while char_window[-1] == "-":
            offset += 1
            char_window = consensus[i : i + offset]
        if char_window.replace("-", "") == "CG":
            # remove this scoring
            # seq1=subfam, seq2=chrom
            # sub_matrix[seq1_char + seq2_char]
            print(char_window)
            print(target[i : i + offset])
        i += offset - 1
    return a.score - remove_score
