from typing import Dict, TextIO

from .constants import PROB_SKIP, PROB_SKIP_TR
from .exceptions import FileFormatException


def read_prior_counts(
    prior_counts_file: TextIO,
    use_trs: bool,
) -> Dict[str, float]:
    subfam_counts: Dict[str, float] = {"skip": PROB_SKIP}

    if use_trs:
        subfam_counts["Tandem Repeat"] = PROB_SKIP_TR
        prob_skip = PROB_SKIP + PROB_SKIP_TR
    else:
        prob_skip = PROB_SKIP

    # Burn the first line headers
    try:
        next(prior_counts_file)
    except StopIteration:
        raise FileFormatException(prior_counts_file.name, 1)

    total_count: float = 0
    current_line_number = 2

    for line in prior_counts_file:
        try:
            subfam, count = line.strip().split()
            subfam_counts[subfam] = float(count)
            total_count += float(count)
        except ValueError:
            raise FileFormatException(
                prior_counts_file.name,
                current_line_number,
            )

        current_line_number += 1

    for key in subfam_counts.keys():
        subfam_counts[key] = (1 - prob_skip) * subfam_counts[key] / total_count

    return subfam_counts
