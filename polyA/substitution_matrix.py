from typing import Dict, TextIO, Tuple

"""
TODO: Explain what this is...
"""
CharacterPositions = Dict[str, int]

"""
TODO: Explain what this is...
"""
SubstitutionMatrix = Dict[int, float]


def load_substitution_matrix(
    file: TextIO,
) -> Tuple[CharacterPositions, SubstitutionMatrix]:
    """
    Load the substitution matrix along with the map of character positions.

    TODO: Explain this better, like what it's used for
    """
    characterPositions: CharacterPositions = {}

    firstLine = next(file)
    for index, character in enumerate(firstLine):
        # TODO: Remove the junk first cuz it messes up the index
        if character == " " or character == "\t":
            continue
        characterPositions[character] = index

    subMatrixColumnCount = len(characterPositions)

    substitutionMatrix: SubstitutionMatrix = {}

    count = 0
    for index, line in enumerate(file):
        subScores = line.split(" ")
        for columnIndex, subScore in enumerate(subScores):
            substitutionMatrix[count * subMatrixColumnCount + columnIndex] = float(
                subScore
            )
        count += 1

    return characterPositions, substitutionMatrix

    # fill that shit
