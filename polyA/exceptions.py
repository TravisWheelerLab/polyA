class FileFormatException(Exception):
    """
    An exception that should be raised whenever we attempt to parse
    an input file and find that it is in the incorrect format.
    """

    path: str
    line_number: int

    def __init__(self, path: str, line_number: int):
        self.path = path
        self._line_number = line_number
