from os import path, mkdir
from sys import stdout
from typing import TextIO, Tuple


class Output:
    """
    A class to manage creating output file objects.
    """

    __base_filename: str
    __output_to_file: bool

    def __init__(self, output_path: str, output_to_file: bool):
        base_filename: str
        if path.isdir(output_path):
            base_filename = path.join(output_path, "output")
        else:
            base_filename = output_path

        if base_filename.endswith(path.sep):
            # We got a directory path but it's not an actual
            # directory so we create it then fix the base
            mkdir(base_filename)
            base_filename = path.join(base_filename, "output")

        self.__base_filename = base_filename
        self.__output_to_file = output_to_file

    def get_heatmap(self, index: int) -> TextIO:
        return open(f"{self.__base_filename}.{index}.heatmap", "w")

    def get_results(self) -> TextIO:
        if not self.__output_to_file:
            return stdout

        results_file = open(f"{self.__base_filename}.results", "w")
        return results_file

    def get_soda(self, index: int) -> Tuple[TextIO, TextIO]:
        viz_file = open(f"{self.__base_filename}.{index}.html", "w")

        return viz_file
