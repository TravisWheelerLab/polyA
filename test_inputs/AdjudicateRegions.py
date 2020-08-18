from getopt import getopt
from math import inf, log
import re
from sys import argv, stdout, stderr
from typing import Dict, List, Tuple
import os
import json

from polyA.calculate_score import calculate_score
from polyA.collapse_matrices import collapse_matrices
from polyA.confidence_cm import confidence_cm
from polyA.fill_align_matrix import fill_align_matrix
from polyA.fill_confidence_matrix import fill_confidence_matrix
from polyA.fill_consensus_position_matrix import fill_consensus_position_matrix
from polyA.fill_probability_matrix import fill_probability_matrix
from polyA.fill_support_matrix import fill_support_matrix
from polyA.load_alignments import load_alignments
from polyA.pad_sequences import pad_sequences
from polyA.printers import print_matrix_support



# -----------------------------------------------------------------------------------#
#			FUNCTIONS													   			#
# -----------------------------------------------------------------------------------#


def GetPath(temp_id: int, columns: List[int], ids: List[int], changes_orig: List[str],
            changes_position_orig: List[int], columns_orig: List[int], subfams_collapse: List[str],
            last_column: List[float], active_cells_collapse: Dict[int, List[int]],
            origin_matrix: Dict[Tuple[int, int], int], same_subfam_change_matrix: Dict[Tuple[int, int], int]) -> Tuple[
    int, List[int], List[str]]:
    """
    using origin matrix, back traces through the 2D array to get the subfam path (most probable
    path through the DP matrix)
    finds where the path switches to a different row and populates Changes and ChangesPosition
    reverses Changes and ChangesPosition because it's a backtrace so they are initially backwards
    jumps over removed/empty columns when necessary

    assigns IDs to each col in matrices (corresponds to a nucleotide position in target/chrom
    sequence) - cols with same ID are part of same initial subfam

    input:
    temp_id: current id number being used - makes it so new ids are unique
    columns: list of non empty columns in matrices
    ids: list of ids for each column in matrices
    changes_orig: original changes from first DP trace
    changes_pos_orig: original changes_position from first DP trace
    subfams_collapse: subfamily names for the rows
    last_column: last column in probability matrix. Used to find where to start backtrace.
    active_cells_collapse: holds which rows haves values for each column
    origin_matrix: origin matrix
    same_subfam_change_matrix: parallel to origin_matrix, if 1 - came from the same subfam, but
    got a change transition probability

    output:
    temp_id: updated current id number being used after function completes
    changes_position: which columns (positions in target/chrom seq) switch to different subfam
    changes: parallel array to changes_position - what subfam is being switches to
    updates input list ids

    >>> non_cols = [0, 1, 2, 3]
    >>> idss = [0, 0, 0, 0]
    >>> subs = ["s1", "s2"]
    >>> active_col = {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]}
    >>> last_col = [-100, -10]
    >>> orig_mat = {(0, 0): 0, (1, 0): 1, (0, 1): 0, (1, 1): 0, (0, 2): 0, (1, 2): 0, (0, 3): 0, (1, 3): 1}
    >>> same_sub_mat = {}
    >>> (temp_idd, changes_pos, changess) = GetPath(1111, non_cols, idss, [], [], [], subs, last_col, active_col, orig_mat, same_sub_mat)
    >>> temp_idd
    3579
    >>> changes_pos
    [0, 2, 4]
    >>> changess
    ['s1', 's2']
    >>> idss
    [2345, 2345, 1111, 1111]
    """

    maxxx: float = -inf
    max_row_index: int = 0

    changes_position: List[int] = []
    changes: List[str] = []

    #which row to start backtrace
    for i in range(len(active_cells_collapse[columns[-1]])):
        if maxxx < last_column[i]:
            maxxx = last_column[i]
            max_row_index = active_cells_collapse[columns[-1]][i]

    prev_row_index: int = origin_matrix[max_row_index, columns[-1]]

    ids[columns[-1]] = temp_id

    changes_position.append(len(columns))

    # already added the last col, but this adds the one before cols so still start at last col
    for columns_index in range(len(columns) - 1, 1, -1):

        prev_column: int = columns[columns_index - 1]
        curr_column: int = columns[columns_index]

        ids[columns[columns_index - 1]] = temp_id

        # updates the original node labels if they change when being stitched
        for i in range(len(changes_position_orig) - 1):
            if columns_orig[changes_position_orig[i]] == prev_column:
                changes_orig[i] = subfams_collapse[origin_matrix[prev_row_index, curr_column]]

        if prev_row_index != origin_matrix[prev_row_index, prev_column]:
            temp_id += 1234
            changes_position.append(columns_index - 1)
            changes.append(subfams_collapse[prev_row_index])
        else:
            if (prev_row_index, prev_column) in same_subfam_change_matrix:
                temp_id += 1234
                changes_position.append(columns_index - 1)
                changes.append(subfams_collapse[prev_row_index])

        prev_row_index = origin_matrix[prev_row_index, prev_column]

    ids[columns[0]] = temp_id
    changes_position.append(0)
    changes.append(subfams_collapse[prev_row_index])

    changes.reverse()
    changes_position.reverse()

    # changes ID for next round of stitching, so when starts stitching will have unique ID
    temp_id += 1234

    return (temp_id, changes_position, changes)


def PrintChanges(columns: List[int], changes: List[str], changes_position: List[int]) -> None:
    """
    just for debugging
    prints out the changes positions and subfams in the matrices, use for each iteration of the node extraction
    """
    i: int = 0
    while i < len(changes):
        stdout.write(str(columns[changes_position[i]]))
        stdout.write("\t")
        stdout.write(str(IDs[columns[changes_position[i]]]))
        stdout.write("\t")
        stdout.write(f"{changes[i]}\n")
        i = i + 1


def PrintNodeConfidence(nodes: int, changes: List[str], subfams_collapse: Dict[str, int],
                        node_confidence: Dict[Tuple[str, int], float]) -> None:
    """
    used for debugging
    prints out matrix that holds node confidence
    """
    for i in range(nodes):
        stdout.write(f"{changes[i]} ")
    stdout.write("\n")

    for subfam in subfams_collapse:
        stdout.write(f"{subfam} ")
        for j in range(nodes):
            if (subfam, j) in node_confidence:
                stdout.write(f"{node_confidence[subfam, j]} ")
            else:
                stdout.write(f"-inf ")
        stdout.write("\n")


def FillNodeConfidence(nodes: int, start_all: int, gap_init: int, gap_ext: int, lamb: float, infilee: str, columns: List[int],
                       starts: List[int], stops: List[int], changes_position: List[int], subfams: List[str], subfam_seqs: List[str], chrom_seqs: List[str],
                       subfam_countss: Dict[str, float], sub_matrix: Dict[str, int]) -> Dict[Tuple[str, int], float]:
    """
    Using changes_position indentified boundaries for all nodes, and computes confidence
    values for each node. First fills matrix with node alignment scores, then reuses matrix
    for confidence scores.

    input:
    all input needed for CalcScore() and confidence_cm()
    nodes: number of nodes
    changes_pos: node boundaries
    columns: all non empty columns in matrices

    output:
    node_confidence: Hash implementation of sparse 2D matrix that holds confidence values
    for whole nodes. Used during stitching process. Tuple[str, int] is key that maps a subfamily
    and node number to a confidence score.

    >>> non_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> strts = [0, 0, 0]
    >>> stps = [10, 10, 10]
    >>> change_pos = [0, 3, 7, 10]
    >>> names = ["skip", "n1", "n2"]
    >>> s_seqs = ['', 'AAA-TTTTT-', 'TTTTTTTTTT']
    >>> c_seqs = ['', 'TTTTTTTTTT', 'TTTTTTTTTT']
    >>> counts = {"skip": .33, "n1": .33, "n2": .33}
    >>> sub_mat = {"AA":1, "AT":-1, "TA":-1, "TT":1}
    >>> node_conf = FillNodeConfidence(3, 0, -25, -5, 0.1227, "infile", non_cols, strts, stps, change_pos, names, s_seqs, c_seqs, counts, sub_mat)
    >>> node_conf
    {('skip', 0): 0.0, ('n1', 0): 0.3751243838973974, ('n2', 0): 0.6248756161026026, ('skip', 1): 0.0, ('n1', 1): 0.09874227070127324, ('n2', 1): 0.9012577292987267, ('skip', 2): 0.0, ('n1', 2): 0.09874227070127327, ('n2', 2): 0.9012577292987267}
    """

    node_confidence_temp: List[float] = [0.0 for _ in range(len(subfams) * nodes)]
    node_confidence: Dict[Tuple[str, int], float] = {}

    #holds all chrom seqs and the offset needed to get to the correct position in the chrom
    #seq - remember gaps in chrom seq are skipped over for matrix position
    chrom_seq_offset: Dict[int, int] = {}

    #first node
    begin_node0: int = columns[changes_position[0]]

    for subfam_index0 in range(1, len(subfams)):

        count: int = 0
        for i in range(begin_node0, columns[changes_position[1]]-columns[changes_position[0]]):
            if chrom_seqs[subfam_index0][i] == '-':
                count += 1

        chrom_seq_offset[subfam_index0] = count

        end_node0: int = columns[changes_position[1]] + count

        subfam0: str = subfam_seqs[subfam_index0][begin_node0:end_node0]
        chrom0: str = chrom_seqs[subfam_index0][begin_node0:end_node0]

        align_score0: float = 0.0
        #if whole alignment is padding - don't run CalcScore
        if end_node0 >= starts[subfam_index0]-start_all and begin_node0 <= stops[subfam_index0] - start_all:
            align_score0 = calculate_score(gap_ext, gap_init, subfam0, chrom0, '', '', sub_matrix)
        node_confidence_temp[subfam_index0 * nodes + 0] = align_score0

    #middle nodes
    for node_index in range(1, nodes-1):

        for subfam_index in range(1, len(subfams)):

            begin_node: int = columns[changes_position[node_index]] + chrom_seq_offset[subfam_index]

            count: int = 0
            for i in range(begin_node, begin_node + columns[changes_position[node_index + 1]] - columns[changes_position[node_index]]):
                if chrom_seqs[subfam_index][i] == '-':
                    count += 1

            chrom_seq_offset[subfam_index] = chrom_seq_offset[subfam_index] + count

            end_node: int = columns[changes_position[node_index+1]] + chrom_seq_offset[subfam_index] + count

            lastprev_subfam: str = subfam_seqs[subfam_index][begin_node-1]
            lastprev_chrom: str = chrom_seqs[subfam_index][begin_node - 1]
            subfam: str = subfam_seqs[subfam_index][begin_node:end_node]
            chrom: str = chrom_seqs[subfam_index][begin_node:end_node]

            align_score: float = 0.0
            # if whole alignment is padding - don't run CalcScore
            if end_node >= starts[subfam_index] - start_all and begin_node <= stops[subfam_index] - start_all:
                align_score = calculate_score(gap_ext, gap_init, subfam, chrom, lastprev_subfam, lastprev_chrom, sub_matrix)

            node_confidence_temp[subfam_index * nodes + node_index] = align_score


    # does last node
    for subfam_index2 in range(1, len(subfams)):
        begin_node2: int = columns[changes_position[-2]] + chrom_seq_offset[subfam_index2]

        count: int = 0
        for i in range(begin_node2,
                       begin_node2 + columns[changes_position[-1] - 1] - columns[changes_position[-2]]):
            if chrom_seqs[subfam_index2][i] == '-':
                count += 1

        chrom_seq_offset[subfam_index2] = chrom_seq_offset[subfam_index2] + count

        end_node2: int = columns[changes_position[-1] - 1] + chrom_seq_offset[subfam_index2] + count

        lastprev_subfam2: str = subfam_seqs[subfam_index2][begin_node2 - 1]
        lastprev_chrom2: str = chrom_seqs[subfam_index2][begin_node2 - 1]

        subfam2: str = subfam_seqs[subfam_index2][begin_node2:end_node2+1]
        chrom2: str = chrom_seqs[subfam_index2][begin_node2:end_node2+1]

        align_score2: float = 0.0
        # if whole alignment is padding - don't run CalcScore
        if end_node2 >= starts[subfam_index2] - start_all and begin_node0 <= stops[subfam_index2] - start_all:
            align_score2 = calculate_score(gap_ext, gap_init, subfam2, chrom2, lastprev_subfam2, lastprev_chrom2,
                                                sub_matrix)

        node_confidence_temp[subfam_index2 * nodes + nodes - 1] = align_score2

    # reuse same matrix and compute confidence scores for the nodes
    for node_index4 in range(nodes):
        temp: List[float] = []
        for row_index in range(1, len(subfams)):
            temp.append(node_confidence_temp[row_index * nodes + node_index4])

        confidence_temp: List[float] = confidence_cm(lamb, infilee, temp, subfam_countss, subfams)

        for row_index2 in range(len(confidence_temp)):
            node_confidence_temp[(row_index2+1) * nodes + node_index4] = confidence_temp[row_index2]

    # collapse node_confidence down same way supportmatrix is collapsed - all seqs of
    # the same subfam are put in the same row
    # not a sparse hash - holds the 0s
    for node_index5 in range(nodes):
        for row_index3 in range(len(subfams)):
            if (subfams[row_index3], node_index5) in node_confidence:
                node_confidence[subfams[row_index3], node_index5] += node_confidence_temp[
                    row_index3 * nodes + node_index5]
            else:
                node_confidence[subfams[row_index3], node_index5] = node_confidence_temp[
                    row_index3 * nodes + node_index5]

    return node_confidence


def PrintPathGraph(nodes: int, changes: List[str], path_graph: List[int]) -> None:
    """
    used for debugging
    prints out the matrix that holds the path graph - 1 means there is an edge, 0 means no edge
    """
    stdout.write(" ")
    for i in range(nodes):
        stdout.write(f"{changes[i]} ")
    stdout.write("\n")

    for i in range(nodes):
        for j in range(nodes):
            stdout.write(f"{path_graph[i * nodes + j]}\t")
        stdout.write(f"{changes[i]}\n")
    stdout.write("\n")


def FillPathGraph(nodes: int, columns: List[int], changes: List[str], changes_position: List[int], consensus_matrix_collapse: Dict[Tuple[int, int], int],
                  strand_matrix_collapse: Dict[Tuple[int, int], str], node_confidence: Dict[Tuple[str, int], float], subfams_collapse_index: Dict[str, int]) -> \
        List[int]:
    """
    finds alternative paths through the nodes - used for stitching to find nested elements
    and to stitch back together elements that have been inserted into.

    Alternative edges only added when confidence of other node label is above a certain
    threshold and when the alignments are relatively contiguous in the consensus sequence.

    input:
    nodes: number of nodes
    columns: list with all non empty columns in matrices
    changes: list of node boundaries
    changes_position: list of subfams identified for each node
    consensus_matrix_collapse: matrix holding alignment position in consensus/subfam seqs -
    used to identify whether nodes are spacially close enough for an alternative edge
    strand_matrix_collapse: matrix holding strand for alignments - nodes can only be connected
    with alternative edge if on same strand
    node_confidence: competing annoation confidence for nodes - nodes can only be connected
    with alternative edge if subfam identified for sink node is a competing annotation in
    source node and is above a certain confidence
    subfams_collapse_index: maps subfam names to their row indices in collapsed matrices

    output:
    path_graph: Flat array reprepsentation of 2D matrix. Graph used during stitching. Maps nodes to all other nodes and holds
    values for if there is an alternative edge between the nodes.

    >>> non_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> changess = ["s1", "s2", "s1"]
    >>> changes_pos = [0, 3, 6, 10]
    >>> con_mat_col = {(0, 0): 0, (0, 1): 1, (0, 2): 2, (0, 3): 3, (0, 6): 6, (0, 7): 7, (0, 8): 8, (0, 9): 9, (1, 3): 1, (1, 4): 2, (1, 5): 3, (1, 6): 4}
    >>> strand_mat_col = {(0, 0): '+', (0, 1): '+', (0, 2): '+', (0, 3): '+', (0, 6): '+', (0, 7): '+', (0, 8): '+', (0, 9): '+', (1, 3): '-', (1, 4): '-', (1, 5): '-', (1, 6): '-'}
    >>> node_conf = {('s1', 0): 0.9, ('s1', 1): 0.5, ('s1', 2): 0.9, ('s1', 0): 0.1, ('s2', 1): 0.5, ('s2', 2): 0.1}
    >>> sub_col_ind = {'s1': 0, 's2': 1}
    >>> FillPathGraph(3, non_cols, changess, changes_pos, con_mat_col, strand_mat_col, node_conf, sub_col_ind)
    [0, 1, 1, 0, 0, 1, 0, 0, 0]
    """

    path_graph: List[int] = []

    for i in range(nodes * nodes):
        path_graph.append(0)

    # filling beginning path graph with straight line through the nodes
    for i in range(nodes - 1):
        path_graph[i * nodes + i + 1] = 1

    for sink_node_index in range(2, nodes):  # don't need to add edges between first 2 nodes because already there
        sink_subfam: str = str(changes[sink_node_index])
        sink_subfam_start: int = consensus_matrix_collapse[subfams_collapse_index[sink_subfam], columns[changes_position[sink_node_index]]]
        sink_strand: str = strand_matrix_collapse[subfams_collapse_index[sink_subfam], columns[changes_position[sink_node_index]]]

        if sink_subfam != "skip":  # don't want to add alternative edges to skip nodes
            for source_node_index in range(sink_node_index - 1):
                source_subfam: str = changes[source_node_index]
                sourceConf: float = node_confidence[sink_subfam, source_node_index]  # sink subfam confidence in source node
                sinkConf: float = node_confidence[source_subfam, sink_node_index]  # source subfam confidence in sink node
                source_subfam_index = subfams_collapse_index[source_subfam]
                source_col: int = columns[changes_position[source_node_index + 1] - 1]

                if (source_subfam_index, source_col) in consensus_matrix_collapse:

                    source_subfam_stop = consensus_matrix_collapse[source_subfam_index, source_col]
                    source_strand = strand_matrix_collapse[source_subfam_index, source_col]

                    # adds in edge if the subfam of the sink is at the source node and if it's
                    # confidence >= 20%, and if the source is before the sink in the consensus sequence

                    # FIXME - not sure what this confidence threshold should be
                    if sourceConf >= 0.2 or sinkConf >= 0.2:
                        if sink_strand == '+' and sink_strand == source_strand:
                            # FIXME- not sure what this overlap should be .. just allowed 50 for now
                            if source_subfam_stop <= sink_subfam_start + 50:
                                path_graph[source_node_index * nodes + sink_node_index] = 1
                        elif sink_strand == '-' and sink_strand == source_strand:
                            if source_subfam_stop + 50 >= sink_subfam_start:
                                path_graph[source_node_index * nodes + sink_node_index] = 1

    return path_graph


def ExtractNodes(num_col: int, nodes: int, columns: List[int], changes_position: List[int],
                 path_graph: List[int]) -> int:
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
    updates columns (NonEmptyColumns)
    updated_num_col: updated number of columns in matrices after nodes have been removed

    >>> non_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> change_pos = [0, 3, 7]
    >>> path_graph = [0, 1, 1, 0, 0, 1, 0, 0, 0]
    >>> ExtractNodes(10, 3, non_cols, change_pos, path_graph)
    10
    >>> non_cols
    [0, 1, 2, 7, 8, 9]
    """

    remove_starts: List[int] = []
    remove_stops: List[int] = []

    updated_num_col: int = num_col

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

    # when removing from the end, have to update num_cols because don't want the end of the matrix anymore
    col_index: int = nodes - 1
    while remove_nodes[col_index]:
        updated_num_col = columns[changes_position[col_index] - 1]
        col_index -= 1

    # removing inserted elements from NonEmptyColumns so they can be ignored
    total: int = 0
    for i in range(len(remove_stops)):
        del columns[remove_starts[i] - total:remove_stops[i] - total]
        # 	helps with offset, when first part is spliced out need an offset to know where to splice out for second part
        total += (remove_stops[i] - remove_starts[i])

    return updated_num_col


def PrintResults(changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int], ids: List[int]) -> None:
    """
    prints the final results
    prints start and stop in terms of position in matrix
    """
    stdout.write("start\tstop\tID\tname\n")
    stdout.write("----------------------------------------\n")
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != 'skip':
            stdout.write(str(columns_orig[changespos_orig[i]]))
            stdout.write("\t")
            stdout.write(str(columns_orig[changespos_orig[i + 1] - 1]))
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def PrintResultsSequence(edgestart: int, changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int],
                         ids: List[int]) -> None:
    """
    prints final results
    prints start and stop in terms of input chrom sequence
    """
    # stdout.write("start\tstop\tID\tname\n")
    # stdout.write("----------------------------------------\n")
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != 'skip':
            stdout.write(str(columns_orig[changespos_orig[i]] + edgestart))
            stdout.write("\t")
            stdout.write(str(columns_orig[changespos_orig[i + 1] - 1] + edgestart))
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def PrintResultsChrom(edgestart: int, chrom_start: int, changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int],
                         ids: List[int]) -> None:
    """
    prints final results
    prints start and stop in terms of chromosome/target position
    """
    stdout.write("start\tstop\tID\tname\n")
    stdout.write("----------------------------------------\n")
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != 'skip':
            stdout.write(str(columns_orig[changespos_orig[i]] + edgestart + chrom_start))
            stdout.write("\t")
            stdout.write(str(columns_orig[changespos_orig[i + 1] - 1] + edgestart + chrom_start))
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def PrintResultsViz(start_all: int, outfile: str, outfile_json: str, chrom: str, chrom_start: int, subfams: List[str], changes_orig: List[str], changes_position_orig: List[int],
                    columns_orig: List[int], consensus_lengths: Dict[str, int],
                    strand_matrix_collapse: Dict[Tuple[int, int], str],
                    consensus_matrix_collapse: Dict[Tuple[int, int], int], subfams_collapse_index: Dict[str, int], node_confidence_orig: Dict[Tuple[str, int], float]):
    """
    prints the results in the proper format to input into the SODA visualization tool

    format described here:
    https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=joinedRmsk&hgta_table=rmskJoinedCurrent&hgta_doSchema=describe+table+schema

    """

    id: int = 0

    length = len(changes_position_orig) - 1

    used: List[int] = [1] * length  # wont print out the results of the same thing twice

    json_dict_id: Dict[str, Dict[str, Dict[str, float]]] = {}

    with open(outfile, 'w') as out:

        i = 0
        while i < length:

            sub_id: int = 0
            json_dict_subid: Dict[str, Dict[str, float]] = {}

            if changes_orig[i] != 'skip' and used[i]:
                subfam: str = changes_orig[i]
                strand: str = strand_matrix_collapse[subfams_collapse_index[subfam], columns_orig[changes_position_orig[i]]]
                consensus_start: int
                consensus_stop: int

                left_flank: int
                right_flank: int

                if strand == "-":
                    left_flank = consensus_lengths[subfam] - consensus_matrix_collapse[
                        subfams_collapse_index[subfam], columns_orig[changes_position_orig[i]]]
                    right_flank = consensus_matrix_collapse[
                                      subfams_collapse_index[subfam], columns_orig[changes_position_orig[i + 1] - 1]] - 1
                else:
                    left_flank = consensus_matrix_collapse[subfams_collapse_index[subfam], columns_orig[changes_position_orig[i]]] - 1
                    right_flank = consensus_lengths[subfam] - consensus_matrix_collapse[
                        subfams_collapse_index[subfam], columns_orig[changes_position_orig[i + 1] - 1]]

                align_start: int = chrom_start + (columns_orig[changes_position_orig[i]] + start_all)
                feature_start: int = align_start - left_flank
                align_stop: int = chrom_start + (columns_orig[changes_position_orig[i + 1] - 1] + start_all)
                feature_stop: int = align_stop + right_flank

                block_start_matrix: int = columns_orig[changes_position_orig[i]]

                block_count: int = 3
                id += 1

                json_dict_subfam_i: Dict[str, float] = {}

                for subfam_i in range(1, len(subfams)):
                    subfamm = subfams[subfam_i]
                    if NodeConfidenceOrig[subfamm, i] > 0.001:
                        json_dict_subfam_i[subfamm] = NodeConfidenceOrig[subfamm, i]

                json_dict_subid[str(id) + "-" + str(sub_id)] = sorted(json_dict_subfam_i.items(), key=lambda x: x[1], reverse=True)
                sub_id += 1

                block_start: List[str] = []
                block_size: List[str] = []

                block_start.append("-1")
                block_start.append(str(left_flank + 1))
                block_start.append("-1")

                block_size.append(str(left_flank))
                block_size.append(
                    str(columns_orig[changes_position_orig[i + 1] - 1] - columns_orig[changes_position_orig[i]] + 1))
                block_size.append(str(right_flank))

                j: int = i + 1
                while j < length:
                    if changes_orig[j] != 'skip':

                        if IDs[columns_orig[changes_position_orig[i]]] == IDs[columns_orig[changes_position_orig[j]]]:

                            if strand == "-":
                                right_flank = consensus_matrix_collapse[
                                                  subfams_collapse_index[changes_orig[j]], columns_orig[changes_position_orig[j + 1] - 1]] - 1
                            else:
                                right_flank = consensus_lengths[subfam] - consensus_matrix_collapse[
                                    subfams_collapse_index[subfam], columns_orig[changes_position_orig[j + 1] - 1]]

                            del block_size[-1]
                            block_size.append("0")

                            align_stop: int = chrom_start + (columns_orig[changes_position_orig[j + 1] - 1] + start_all)
                            feature_stop: int = align_stop + right_flank

                            block_start.append(
                                str(columns_orig[changes_position_orig[j]] + 1 - block_start_matrix + left_flank))
                            block_start.append("-1")

                            block_size.append(str(columns_orig[changes_position_orig[j + 1] - 1] - columns_orig[
                                changes_position_orig[j]] + 1))

                            block_size.append(str(right_flank))

                            block_count += 2

                            used[j] = 0

                            json_dict_subfam_j: Dict[str, float] = {}

                            for subfam_i in range(1, len(Subfams)):
                                subfamm = Subfams[subfam_i]
                                if NodeConfidenceOrig[subfamm, j] > 0.001:
                                    json_dict_subfam_j[subfamm] = node_confidence_orig[subfamm, j]

                            json_dict_subid[str(id) + "-" + str(sub_id)] = sorted(json_dict_subfam_j.items(), key=lambda x: x[1], reverse=True)
                            sub_id += 1

                    j += 1

                json_dict_id[str(id)] = json_dict_subid

                out.write("000 " + chrom + " " + str(feature_start) + " " + str(
                    feature_stop) + " " + subfam + " 0 " + strand + " " + str(align_start) + " " + str(
                    align_stop) + " 0 " + str(block_count) + " " + (",".join(block_size)) + " " + (
                              ",".join(block_start)) + " " + str(id))
                out.write("\n")

            used[i] = 0
            i += 1

    #prints json file with confidence values for each annotation
    with open(outfile_json, 'w') as out_json:
        out_json.write(json.dumps(json_dict_id))



# -----------------------------------------------------------------------------------#
#			   MAIN														   			#
# -----------------------------------------------------------------------------------#


if __name__ == "__main__":

    GapInit: int = -25
    GapExt: int = -5
    Lamb: float = 0.0
    EslPath = ""
    ChunkSize: int = 31
    SameProbLog: float = 0.0
    ChangeProb: float = 10 ** -45
    ChangeProbLog: float = 0.0  # Reassigned later
    ChangeProbSkip: float = 0.0  # Reassigned later
    SameProbSkip: float = 0.0
    SkipAlignScore: float = 0.0

    StartAll: int = 0  # Reassigned later
    StopAll: int = 0  # Reassigned later
    ID: int = 1111

    infile_prior_counts: str = ""
    outfile_viz: str = ""
    outfile_heatmap: str = ""

    help: bool = False  # Reassigned later
    prin: bool = False  # Reassigned later
    printMatrixPos: bool = False  # Reassigned later
    printSeqPos: bool = False

    helpMessage: str = f"""
    usage: {argv[0]} alignFile subMatrixFile\n
    ARGUMENTS
        --GapInit[-25]
        --getExt[-5]
        --lambda [will calc from matrix if not included]
        --eslPath [specify path to easel]
        --segmentsize (must be odd) [31]
        --changeprob[1e-45]
        --priorCounts PriorCountsFile
    
    OPTIONS
        --help - display help message
        --matrixpos - prints subfam changes in matrix position instead of genomic position
        --sequencepos - prints subfam changes in sequence position instead of genomic position
        --viz outfile - prints output format for SODA visualization
        --heatmap outfile - prints probability file for input into heatmap
    """

    raw_opts, args = getopt(argv[1:], "", [
        "GapInit=",
        "GapExt=",
        "skipScore=",
        "lambda=",
        "eslPath=",
        "segmentsize=",
        "changeprob=",
        "priorCounts=",
        "viz=",
        "heatmap=",

        "help",
        "matrixpos",
        "seqpos",
    ])
    opts = dict(raw_opts)

    GapInit = int(opts["--GapInit"]) if "--GapInit" in opts else GapInit
    GapExt = int(opts["--GapExt"]) if "--GapExt" in opts else GapExt
    Lamb = float(opts["--lambda"]) if "--lambda" in opts else Lamb
    EslPath = str(opts["--eslPath"]) if "--eslPath" in opts else EslPath
    ChunkSize = int(opts["--segmentsize"]) if "--segmentsize" in opts else ChunkSize
    ChangeProb = float(opts["--changeprob"]) if "--changeprob" in opts else ChangeProb
    infile_prior_counts = str(opts["--priorCounts"]) if "--priorCounts" in opts else infile_prior_counts
    outfile_viz = str(opts["--viz"]) if "--viz" in opts else outfile_viz
    outfile_heatmap = str(opts["--heatmap"]) if "--heatmap" in opts else outfile_heatmap
    help = "--help" in opts
    printMatrixPos = "--matrixpos" in opts
    printSeqPos = "--seqpos" in opts

    outfile_conf = outfile_viz + ".json"

    if help:
        print(helpMessage)
        exit(0)

    # input is alignment file of hits region and substitution matrix
    infile: str = args[0]
    infile_matrix: str = args[1]

    # Other open was moved down to where we load the alignments file
    with open(infile_matrix) as _infile_matrix:
        in_matrix: List[str] = _infile_matrix.readlines()

    #if command line option included to use prior counts into in confidence calculations
    if infile_prior_counts:
        with open(infile_prior_counts) as _infile_prior_counts:
            in_counts: List[str] = _infile_prior_counts.readlines()

    #if lambda isn't included at command line, run esl_scorematrix to calculate it from scorematrix
    if not Lamb:
        esl_stream = os.popen(EslPath + 'esl_scorematrix --dna 25p41g_edited.matrix')
        esl_output = esl_stream.read()
        esl_output_list = re.split(r"\n+", esl_output)
        lambda_list = re.split(r"\s+", esl_output_list[1])
        Lamb = float(lambda_list[2])

    # reads in the score matrix from file and stores in dict that maps 'char1char2' to the score from the
    # input substitution matrix - ex: 'AA' = 8
    SubMatrix: Dict[str, int] = {}
    line = in_matrix[0]
    line = re.sub(r"^\s+", "", line)
    line = re.sub(r"\s+$", "", line)
    chars = re.split(r"\s+", line)

    count: int = 0
    for line in in_matrix[1:]:
        line = re.sub(r"^\s+", "", line)
        line = re.sub(r"\s+$", "", line)
        subScores = re.split(r"\s+", line)
        for i in range(len(subScores)):
            SubMatrix[chars[count] + chars[i]] = int(subScores[i])
        count += 1
    SubMatrix['..'] = 0

    # maps subfam names to genomic prior_count/total_in_genome from input file
    # used during confidence calculations
    SubfamCounts: Dict[str, float] = {}
    PriorTotal: float = 0
    prob_skip = 0.4  # about 60% of genome is TE derived
    if infile_prior_counts:
        for line in in_counts[1:]:
            line = re.sub(r"\n", "", line)
            info = re.split(r"\s+", line)
            SubfamCounts[info[0]] = int(info[1])
            PriorTotal += float(info[1])

        for key in SubfamCounts:
            SubfamCounts[key] = (1 - prob_skip) * SubfamCounts[key] / PriorTotal

        SubfamCounts["skip"] = prob_skip

    Subfams: List[str] = []
    Chroms: List[str] = []
    Scores: List[int] = []
    Strands: List[str] = []
    Starts: List[int] = []
    Stops: List[int] = []
    ConsensusStarts: List[int] = []
    ConsensusStops: List[int] = []
    SubfamSeqs: List[str] = []
    ChromSeqs: List[str] = []
    Flanks: List[int] = []

    AlignMatrix: Dict[Tuple[int, int], float] = {}
    SingleAlignMatrix: Dict[Tuple[int, int], int] = {}
    ConfidenceMatrix: Dict[Tuple[int, int], float] = {}
    SupportMatrix: Dict[Tuple[int, int], float] = {}
    OriginMatrix: Dict[Tuple[int, int], int] = {}
    ConsensusMatrix: Dict[Tuple[int, int], int] = {}
    SameSubfamChangeMatrix: Dict[Tuple[int, int], int] = {}
    ProbMatrixLastColumn: List[float] = []

    NonEmptyColumns: List[int] = []
    ActiveCells: Dict[int, List[int]] = {}

    Changes: List[str] = []
    ChangesPosition: List[int] = []

    IDs: List[int] = []
    ChangesOrig: List[str] = []
    ChangesPositionOrig: List[int] = []
    NonEmptyColumnsOrig: List[int] = []

    # for graph/node part
    NumNodes: int = 0
    NodeConfidence: Dict[Tuple[str, int], float] = {}
    NodeConfidenceOrig: Dict[Tuple[str, int], float] = {}
    PathGraph: List[int] = []
    total: int = 0
    loop: int = 1

    # opens alignment file and stores all Subfams, Scores, Starts, Stops, subfams seqs and Chrom seqs in arrays
    numseqs: int = 0
    with open(infile) as _infile:
        alignments = load_alignments(_infile)
        for alignment in alignments:
            numseqs += 1

            Subfams.append(alignment.subfamily)
            Chroms.append(alignment.chrom)
            Scores.append(alignment.score)
            Strands.append(alignment.strand)
            Starts.append(alignment.start)
            Stops.append(alignment.stop)
            ConsensusStarts.append(alignment.consensus_start)
            ConsensusStops.append(alignment.consensus_stop)
            SubfamSeqs.append(alignment.subfamily_sequence)
            ChromSeqs.append(alignment.sequence)
            Flanks.append(alignment.flank)

    # if there is only one subfam in the alignment file, no need to run anything because we know
    # that subfam is the annotation
    # numseqs = 2 because of the skip state
    if numseqs == 2:
        if printMatrixPos:
            stdout.write("start\tstop\tID\tname\n")
            stdout.write("----------------------------------------\n")
            stdout.write(f"{0}\t{Stops[1] - Starts[1]}\t1111\t{Subfams[1]}\n")
        else:
            stdout.write("start\tstop\tID\tname\n")
            stdout.write("----------------------------------------\n")
            stdout.write(f"{Starts[1]}\t{Stops[1]}\t1111\t{Subfams[1]}\n")
        exit()

    match = re.search(r"(.+):(\d+)-(\d+)", Chroms[1])
    Chrom: str = match.groups()[0]
    ChromStart: int = int(match.groups()[1])
    ChromEnd: int = int(match.groups()[2])
    TargetLen: int = ChromEnd - ChromStart

    #bail out if target sq is < 25 nucls
    #warning if less than 1000 nucls
    #warning if no chrom info given - okay for artificial seq inputs
    if TargetLen == 0:
        stderr.write("WARNING - No chromosome position information given.\nThis is okay if running on artificial sequences, but cannot use command line options --viz or --heatmap.\n\n")
    elif TargetLen <= 25:
        stderr.write("ERROR - Target sequence length needs to be > 25 nucleotides.\n\n")
        exit()
    elif TargetLen < 1000:
        stderr.write("WARNING - Did you mean to run this on a target region < 1000 nucleotides?\n\n")

    ChangeProbLog = log(ChangeProb / (numseqs - 1))
    ChangeProbSkip = (ChangeProbLog / 2)  # jumping in and then out of the skip state counts as 1 jump
    SameProbSkip = ChangeProbLog / 20  # 5% of the jump penalty, staying in skip state for 20nt "counts" as one jump

    # precomputes number of rows in matrices
    rows: int = len(Subfams)
    cols: int = 0  # assign cols in FillAlignMatrix

    # precomputes consensus seq length for PrintResultsViz()
    ConsensusLengths: Dict[str, int] = {}
    if outfile_viz:
        for i in range(1, len(Flanks)):
            if Strands[i] == "+":
                ConsensusLengths[Subfams[i]] = ConsensusStops[i] + Flanks[i]
            else:
                ConsensusLengths[Subfams[i]] = ConsensusStarts[i] + Flanks[i]

    StartAll, StopAll = pad_sequences(ChunkSize, Starts, Stops, SubfamSeqs, ChromSeqs)

    (cols, AlignMatrix) = fill_align_matrix(StartAll, ChunkSize, GapExt, GapInit, SkipAlignScore, SubfamSeqs,
                                          ChromSeqs, Starts, SubMatrix)

    (NonEmptyColumns, ActiveCells, ConsensusMatrix) = fill_consensus_position_matrix(cols, rows, StartAll, SubfamSeqs, ChromSeqs, Starts, Stops, ConsensusStarts, Strands)

    ConfidenceMatrix = fill_confidence_matrix(Lamb, infile_prior_counts, NonEmptyColumns, SubfamCounts, Subfams, ActiveCells,
                                            AlignMatrix)

    SupportMatrix = fill_support_matrix(rows, ChunkSize, StartAll, NonEmptyColumns, Starts, Stops, ConfidenceMatrix)

    collapsed_matrices = collapse_matrices(rows, NonEmptyColumns, Subfams, Strands, ActiveCells, SupportMatrix, ConsensusMatrix)

    SupportMatrixCollapse = collapsed_matrices.support_matrix
    SubfamsCollapse = collapsed_matrices.subfamilies
    ActiveCellsCollapse = collapsed_matrices.active_rows
    ConsensusMatrixCollapse = collapsed_matrices.consensus_matrix
    StrandMatrixCollapse = collapsed_matrices.strand_matrix
    SubfamsCollapseIndex = collapsed_matrices.subfamily_indices

    #if command line option included to output support matrix for heatmap
    if outfile_heatmap:
        with open(outfile_heatmap, "w") as outfile:
            print_matrix_support(cols, StartAll, ChromStart, SupportMatrixCollapse, SubfamsCollapse, file = outfile)

    (ProbMatrixLastColumn, OriginMatrix, SameSubfamChangeMatrix) = fill_probability_matrix(SameProbSkip, SameProbLog, ChangeProbLog,
                                                                               ChangeProbSkip,
                                                                               NonEmptyColumns,
                                                                                         collapsed_matrices,)

    #IDs for each nucleotide will be assigned during DP backtrace
    IDs = [0] * cols

    (ID, ChangesPosition, Changes) = GetPath(ID, NonEmptyColumns, IDs, ChangesOrig, ChangesPositionOrig,
                                             NonEmptyColumnsOrig, SubfamsCollapse, ProbMatrixLastColumn, ActiveCellsCollapse,
                                             OriginMatrix, SameSubfamChangeMatrix)

    # keep the original annotation for reporting results
    ChangesOrig = Changes.copy()
    ChangesPositionOrig = ChangesPosition.copy()
    NonEmptyColumnsOrig = NonEmptyColumns.copy()

    # Finding inserted elements and stitching original elements
    # Steps-
    # 1. Find confidence for nodes
    # 2. Create path graph - Find alternative paths through graph and add those edges
    # 3. Extract all nodes (from dp matrix) that have a single incoming and a single outgoing edge
    # 5. Annotate again with removed nodes
    #   ** stop when all nodes have incoming and outgoing edges <= 1 or there are <= 2 nodes left

    count: int = 0
    while (True):
        count += 1
        NumNodes = len(Changes)

        # breakout of loop if there are 2 or less nodes left
        if (NumNodes <= 2):
            break

        NodeConfidence.clear() #reuse old NodeConfidence matrix

        NodeConfidence = FillNodeConfidence(NumNodes, StartAll, GapInit, GapExt, Lamb, infile_prior_counts, NonEmptyColumns, Starts, Stops, ChangesPosition, Subfams, SubfamSeqs, ChromSeqs, SubfamCounts, SubMatrix)

        # store original node confidence for reporting results
        if count == 1:
            NodeConfidenceOrig = NodeConfidence.copy()

        PathGraph.clear() #reuse old PathGraph
        PathGraph = FillPathGraph(NumNodes, NonEmptyColumns, Changes, ChangesPosition,
                                  ConsensusMatrixCollapse, StrandMatrixCollapse, NodeConfidence, SubfamsCollapseIndex)

        # test to see if there are nodes in the graph that have more that one incoming or outgoing edge,
        # if so - keep looping, if not - break out of the loop
        # if they are all 0, break out of the loop
        test: bool = False
        j: int = 0
        while j < NumNodes:
            i: int = 0
            while i < j - 1:
                if (PathGraph[i * NumNodes + j] == 1):
                    test = True
                i += 1
            j += 1

        if not test:
            break

        cols = ExtractNodes(cols, NumNodes, NonEmptyColumns, ChangesPosition, PathGraph)

        # run DP calculations again with nodes corresponding to inserted elements removed
        # ignores removed nodes because they are no longer in NonEmptyColumns
        (ProbMatrixLastColumn, OriginMatrix, SameSubfamChangeMatrix) = fill_probability_matrix(SameProbSkip, SameProbLog,
                                                                                   ChangeProbLog, ChangeProbSkip,
                                                                                   NonEmptyColumns,
                                                                                   ActiveCellsCollapse,
                                                                                               collapsed_matrices,)

        Changes.clear()
        ChangesPosition.clear()

        (ID, ChangesPosition, Changes) = GetPath(ID, NonEmptyColumns, IDs, ChangesOrig, ChangesPositionOrig,
                                                 NonEmptyColumnsOrig, SubfamsCollapse, ProbMatrixLastColumn, ActiveCellsCollapse,
                                                 OriginMatrix, SameSubfamChangeMatrix)


    #prints results
    if printMatrixPos:
        PrintResults(ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs)
    elif printSeqPos:
        PrintResultsSequence(StartAll, ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs)
    else:
        PrintResultsChrom(StartAll, ChromStart, ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs)

    if outfile_viz:
        PrintResultsViz(StartAll, outfile_viz, outfile_conf, Chrom, ChromStart, Subfams, ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, ConsensusLengths,
                        StrandMatrixCollapse, ConsensusMatrixCollapse, SubfamsCollapseIndex, NodeConfidenceOrig)

