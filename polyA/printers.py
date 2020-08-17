from math import inf
from sys import stdout
from typing import Dict, List, TextIO, Tuple, Union

from polyA.matrices import SupportMatrix


def print_matrix_hash(
    num_col: int,
    num_row: int,
    subfams: List[str],
    matrix: Dict[Tuple[int, int], Union[float, int, str]],
    file: TextIO = stdout,
) -> None:
    """
    Print values contained in a non-collapsed matrix.
    This function exists for debugging purposes only.
    """
    file.write("\t")
    j: int = 0
    while j < num_col:
        file.write(f"{j}\t")
        j += 1
    file.write("\n")

    i: int = 0
    while i < num_row:
        file.write(f"{subfams[i]}\t")
        j: int = 0
        while j < num_col:
            if (i, j) in matrix:
                file.write(f"{matrix[i, j]}")
            else:
                file.write(f"{-inf}")
            file.write("\t")
            j += 1
        file.write("\n")
        i += 1


def print_matrix_support(
    num_col: int,
    start_all: int,
    chrom_start: int,
    matrix: SupportMatrix,
    subfams_collapse: List[str],
    file: TextIO = stdout,
) -> None:
    """
    Prints the given support matrix to `file` in a format appropriate
    to produce a heatmap.
    """
    file.write("\t")

    start: int = chrom_start + start_all
    j: int = 0
    while j < num_col:
        file.write(f"{start}\t")
        start += 1
        j += 1
    file.write("\n")

    for k in range(len(subfams_collapse)):
        file.write(f"{subfams_collapse[k]}\t")
        j: int = 0
        while j < num_col:
            if (k, j) in matrix:
                file.write(str(matrix[k, j]))
            else:
                file.write(f"{-inf}")
            file.write("\t")
            j += 1
        file.write("\n")
