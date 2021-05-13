from os import path, mkdir
from sys import stdout
from typing import TextIO, Tuple


class Output:
    """
    A class to manage creating output file objects.

    TODO(George): Use this for the log file
    """

    __base_filename: str

    def __init__(self, output_path: str):
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

    def get_heatmap(self, index: int) -> TextIO:
        return open(f"{self.__base_filename}.{index}.heatmap", "w")

    def get_results(self) -> TextIO:
        return stdout
        # TODO: For now we use stdout, but later make it an option
        # results_file = open(f"{self.__base_filename}.output", "w")
        # return results_file

    def get_soda(self, index: int) -> Tuple[TextIO, TextIO]:
        viz_file = open(f"{self.__base_filename}.{index}.viz", "w")
        conf_file = open(f"{self.__base_filename}.{index}.viz.json", "w")

        return viz_file, conf_file
