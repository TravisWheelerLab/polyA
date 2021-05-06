from typing import List

from .performance import timeit


@timeit()
def extract_nodes(
    nodes: int,
    columns: List[int],
    changes_position: List[int],
    path_graph: List[int],
) -> None:
    """
    finds nodes that only have one (or less) incoming and one (or less) outgoing edge and adds
    them to remove_starts and remove_stops so they can be extracted from the alignment - all columns
    in removed nodes are removed from NonEmptyColumns (columns) so during next round of DP calculations,
    those nodes are ignored and surrounding nodes can be stitched if necessary,

    input:
    num_col: number of columns in matrices
    nodes: number of nodes in graph
    columns: non empty columns in matrices
    changes_position: node boundaries
    path_graph: 2D array that represents edges in the graph

    output:
    updated_num_col: updated number of columns after nodes are extracted

    >>> non_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> change_pos = [0, 3, 7]
    >>> path_graph = [0, 1, 1, 0, 0, 1, 0, 0, 0]
    >>> extract_nodes(3, non_cols, change_pos, path_graph)
    >>> non_cols
    [0, 1, 2, 7, 8, 9]
    """

    remove_starts: List[int] = []
    remove_stops: List[int] = []

    # boolean for which nodes will be removed
    remove_nodes: List[bool] = [False for _ in range(nodes)]

    # extracting nodes that only have one incoming and one outgoing edge
    num_edges_in: List[int] = [0 for _ in range(nodes)]
    num_edges_out: List[int] = [0 for _ in range(nodes)]

    for row in range(nodes):
        for col in range(nodes):
            num_edges_in[col] += path_graph[row * nodes + col]
            num_edges_out[row] += path_graph[row * nodes + col]

    for node in range(nodes - 1):
        if num_edges_in[node] <= 1 and num_edges_out[node] <= 1:
            remove_starts.append(changes_position[node])
            remove_stops.append(changes_position[node + 1])
            remove_nodes[node] = True

    # deals with last node, so when NumNodes-1 the last remove stop is the end of the matrix
    if num_edges_in[nodes - 1] <= 1 and num_edges_out[nodes - 1] <= 1:
        remove_starts.append(changes_position[nodes - 1])
        remove_stops.append(len(columns) - 1)
        remove_nodes[nodes - 1] = True

    # removing inserted elements from NonEmptyColumns so they can be ignored
    total: int = 0
    for i in range(len(remove_stops)):
        del columns[remove_starts[i] - total : remove_stops[i] - total]
        # 	helps with offset, when first part is spliced out need an offset to know where to splice out for second part
        total += remove_stops[i] - remove_starts[i]
