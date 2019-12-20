from typing import List

from .fill_prob_matrix import fill_prob_matrix
from .serialize_prob_matrix import serialize_prob_matrix
from .support_matrix import deserialize_support_matrix

class Options:
    support_matrix_path: str

    def __init__(self, args: List[str]):
        # TODO: Properly parse the arguments
        self.support_matrix_path = args[0]

def run(options: Options) -> None:
    with open(options.support_matrix_path) as supportMatrixFile:
        supportMatrixLines = supportMatrixFile.readlines()
    
    supportMatrix = deserialize_support_matrix(supportMatrixLines)
    probMatrix, originMatrix = fill_prob_matrix(supportMatrix, list(range(902)), 135)
    serialize_prob_matrix(probMatrix, originMatrix)

