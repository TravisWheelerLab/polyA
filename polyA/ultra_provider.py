import json
import os
from logging import Logger
from typing import Any, Callable, Dict, List, NamedTuple, Tuple

from .performance import timeit


class TandemRepeat(NamedTuple):
    consensus: str
    start: int
    length: int
    stop: int
    position_scores: Tuple[float, ...]

    @staticmethod
    def from_json(json_map: Dict[str, Any]):
        raw_scores = json_map["PositionScoreDelta"].split(":")
        position_scores = []
        for score in raw_scores:
            if abs(float(score)) > 100:
                position_scores.append(0.0)
                Logger(__name__).warning(
                    """
                    Unreasonably-large score detected from ULTRA, 
                    in sequence region {}..{}, score={}
                    """.format(
                        int(json_map["Start"]),
                        int(json_map["Start"]) + int(json_map["Length"]) - 1,
                        score,
                    )
                )
            else:
                position_scores.append(float(score))

        return TandemRepeat(
            consensus=json_map["Consensus"],
            start=int(json_map["Start"]),
            length=int(json_map["Length"]),
            stop=int(json_map["Start"]) + int(json_map["Length"]) - 1,
            position_scores=tuple(position_scores),
        )


class UltraOutput(NamedTuple):
    tandem_repeats: List[TandemRepeat]

    @staticmethod
    def from_json(json_map: Dict[str, Any]):
        repeats = list(map(TandemRepeat.from_json, json_map["Repeats"]))
        return UltraOutput(tandem_repeats=repeats)

    @property
    def tr_count(self) -> int:
        return len(self.tandem_repeats)


UltraProvider = Callable[[], UltraOutput]


class ApplicationUltraProvider:
    _sequence_path: str
    _ultra_output_path: str
    _ultra_path: str

    def __init__(
        self,
        sequence_path: str = "",
        ultra_output_path: str = "",
        ultra_path: str = "ultra",
    ):
        if ultra_output_path == "" and sequence_path == "":
            raise RuntimeError(
                "must provide either ultra output path or sequence file path"
            )

        self._sequence_path = sequence_path
        self._ultra_output_path = ultra_output_path
        self._ultra_path = ultra_path

    @timeit()
    def __call__(self) -> UltraOutput:
        if self._sequence_path:
            ultra_stream = os.popen(
                self._ultra_path + " -ss " + self._sequence_path
            )
            raw_output = json.load(ultra_stream)
        else:
            with open(self._ultra_output_path, "r") as output_file:
                raw_output = json.load(output_file)

        ultra_output = UltraOutput.from_json(raw_output)

        return ultra_output
