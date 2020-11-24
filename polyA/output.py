from os import path
from typing import TextIO, Tuple


class Output:
    """
    A class to manage creating output file objects.
    """
    base_filename: str

    def __init__(self, output_path: str):
        base_filename: str
        if path.isdir(output_path):
            base_filename = path.join(output_path, "output")
        else:
            base_filename = output_path

        if base_filename.endswith(path.sep):
            # We got a directory path but it's not an actual
            # directory so we just quit because we don't know
            # with certainty what the user wanted.
            raise RuntimeError(f"{base_filename} does not exist")

        self.base_filename = base_filename

    def get_heatmap(self, index: int) -> TextIO:
        return open(f"{self.base_filename}.{index}.heatmap", "w")

    def get_soda(self, index: int) -> Tuple[TextIO, TextIO]:
        viz_file = open(f"{self.base_filename}.{index}.viz", "w")
        conf_file = open(f"{self.base_filename}.{index}.viz.json", "w")

        return viz_file, conf_file
