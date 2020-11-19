from typing import Dict, TextIO

from polyA import PROB_SKIP, PROB_SKIP_TR


def read_prior_counts(
    prior_counts_file: TextIO, use_trs: bool
) -> Dict[str, float]:
    subfam_counts: Dict[str, float] = {"skip": PROB_SKIP}

    if use_trs:
        subfam_counts["Tandem Repeat"] = PROB_SKIP_TR
        prob_skip = PROB_SKIP + PROB_SKIP_TR
    else:
        prob_skip = PROB_SKIP

    # TODO: Add some validation for the file

    # Burn the first line headers
    next(prior_counts_file)

    total_count = 0

    for line in prior_counts_file:
        subfam, count = line.strip().split()
        subfam_counts[subfam] = float(count)
        total_count += float(count)

    for key in subfam_counts:
        subfam_counts[key] = (1 - prob_skip) * subfam_counts[key] / total_count

    return subfam_counts
