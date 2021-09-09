class FileFormatException(Exception):
    """
    An exception that should be raised whenever we attempt to parse
    an input file and find that it is in the incorrect format.
    """

    path: str
    line_number: int
    message: str

    def __init__(self, path: str, line_number: int, message: str = ""):
        self.path = path
        self.line_number = line_number
        self.message = message
