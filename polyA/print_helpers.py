from typing import Dict, List

from .alignment import Alignment

ConsensusLengths = Dict[str, int]


def find_consensus_lengths(alignments: List[Alignment]) -> ConsensusLengths:
    """
    Consensus sequence lengths used by the SODA printer.
    """
    consensus_lengths = {}
    for a in alignments:
        flank = a.flank

        if a.strand == "+":
            stop_position = a.consensus_stop
            consensus_lengths[a.subfamily] = stop_position + flank
        else:
            start_position = a.consensus_start
            consensus_lengths[a.subfamily] = start_position + flank

    return consensus_lengths
