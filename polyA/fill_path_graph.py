from typing import Dict, List, Tuple

from polyA.performance import timeit


@timeit()
def fill_path_graph(
    nodes: int,
    columns: List[int],
    changes: List[str],
    changes_position: List[int],
    consensus_matrix_collapse: Dict[Tuple[int, int], int],
    strand_matrix_collapse: Dict[Tuple[int, int], str],
    node_confidence: Dict[Tuple[str, int], float],
    subfams_collapse_index: Dict[str, int],
) -> List[int]:
    """
    finds alternative paths through the nodes - used for stitching to find nested elements
    and to stitch back together elements that have been inserted into.

    Alternative edges only added when on same strand, and confidence of other node label is above a certain
    threshold and when the alignments are relatively contiguous in the consensus sequence.

    input:
    nodes: number of nodes
    columns: list with all non empty columns in matrices
    changes: list of node boundaries
    changes_position: list of subfams identified for each node
    consensus_matrix_collapse: matrix holding alignment position in consensus/subfam seqs -
    used to identify whether nodes are spacially close enough for an alternative edge
    strand_matrix_collapse: matrix holding strand for alignments
    node_confidence: competing annoation confidence for nodes
    subfams_collapse_index: maps subfam names to their row indices in collapsed matrices

    output:
    path_graph: Flat array reprepsentation of 2D matrix. Graph used during stitching. Maps nodes to all other nodes.
    1 if nodes are connected, 0 if not

    >>> non_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> changess = ["s1", "s2", "s1"]
    >>> changes_pos = [0, 3, 6, 10]
    >>> con_mat_col = {(0, 0): 0, (0, 1): 1, (0, 2): 2, (0, 3): 3, (0, 6): 6, (0, 7): 7, (0, 8): 8, (0, 9): 9, (1, 3): 1, (1, 4): 2, (1, 5): 3, (1, 6): 4}
    >>> strand_mat_col = {(0, 0): '+', (0, 1): '+', (0, 2): '+', (0, 3): '+', (0, 6): '+', (0, 7): '+', (0, 8): '+', (0, 9): '+', (1, 3): '-', (1, 4): '-', (1, 5): '-', (1, 6): '-'}
    >>> node_conf = {('s1', 0): 0.9, ('s1', 1): 0.5, ('s1', 2): 0.9, ('s1', 0): 0.1, ('s2', 1): 0.5, ('s2', 2): 0.1}
    >>> sub_col_ind = {'s1': 0, 's2': 1}
    >>> fill_path_graph(3, non_cols, changess, changes_pos, con_mat_col, strand_mat_col, node_conf, sub_col_ind)
    [0, 1, 1, 0, 0, 1, 0, 0, 0]
    """
    path_graph: List[int] = []

    for i in range(nodes * nodes):
        path_graph.append(0)

    # filling beginning path graph with straight line through the nodes
    for i in range(nodes - 1):
        path_graph[i * nodes + i + 1] = 1

    for sink_node_index in range(
        2, nodes
    ):  # don't need to add edges between first 2 nodes because already there
        sink_subfam: str = str(changes[sink_node_index])
        sink_subfam_start: int = consensus_matrix_collapse[
            subfams_collapse_index[sink_subfam],
            columns[changes_position[sink_node_index]],
        ]
        sink_strand: str = strand_matrix_collapse[
            subfams_collapse_index[sink_subfam],
            columns[changes_position[sink_node_index]],
        ]

        if (
            sink_subfam != "skip" and sink_subfam != "Tandem Repeat"
        ):  # don't want to add alternative edges to skip nodes or tandem repeats
            for source_node_index in range(sink_node_index - 1):
                source_subfam: str = changes[source_node_index]
                if source_subfam != "Tandem Repeat":
                    source_conf: float = 0.0
                    sink_conf: float = 0.0
                    if (sink_subfam, source_node_index) in node_confidence:
                        source_conf = node_confidence[
                            sink_subfam, source_node_index
                        ]  # sink subfam confidence in source node
                    if (source_subfam, sink_node_index) in node_confidence:
                        sink_conf = node_confidence[
                            source_subfam, sink_node_index
                        ]  # source subfam confidence in sink node
                    source_subfam_index = subfams_collapse_index[source_subfam]
                    source_col: int = columns[
                        changes_position[source_node_index + 1] - 1
                    ]

                    if (
                        source_subfam_index,
                        source_col,
                    ) in consensus_matrix_collapse:

                        source_subfam_stop = consensus_matrix_collapse[
                            source_subfam_index, source_col
                        ]
                        source_strand = strand_matrix_collapse[
                            source_subfam_index, source_col
                        ]

                        # adds in edge if the subfam of the sink is at the source node and if it's
                        # confidence >= 20%, and if the source is before the sink in the consensus sequence

                        if source_conf >= 0.01 or sink_conf >= 0.01:
                            if (
                                sink_strand == "+"
                                and sink_strand == source_strand
                            ):
                                if source_subfam_stop <= sink_subfam_start + 50:
                                    path_graph[
                                        source_node_index * nodes
                                        + sink_node_index
                                    ] = 1
                            elif (
                                sink_strand == "-"
                                and sink_strand == source_strand
                            ):
                                if source_subfam_stop + 50 >= sink_subfam_start:
                                    path_graph[
                                        source_node_index * nodes
                                        + sink_node_index
                                    ] = 1

    return path_graph
