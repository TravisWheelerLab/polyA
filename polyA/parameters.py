from typing import NamedTuple, TextIO, Union, Callable

EaselRunner = Callable[[], None]


# TODO: Read and write JSON for repeatability
# TODO: Separate the IO stuff into another class
class Parameters(NamedTuple):
    gap_init: int = -25
    gap_ext: int = -5

    input_lambda: float = 0.0
    easel_path: str = ""

    chunk_size: int = 31

    same_prob_log: float = 0.0
    same_prob_skip: float = 0.0

    change_prob: float = 10 ** -45
    change_prob_log: float = 0.0
    change_prob_skip: float = 0.0

    skip_align_score: float = 0.0

    start_all: int = 0
    stop_all: int = 0
    id: int = 1111

    prior_counts_path: str = ""
    soda_path: str = ""
    heatmap_path: str = ""
