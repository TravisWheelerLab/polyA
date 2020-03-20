from typing import Dict, TextIO, Tuple

"""
A dict that associates characters from score matrix file with
positions in score matrix alignment. The keys are characters
and the values as indices into the array.   

TODO: Verify the explanation above and clarify if needed
"""
CharacterPositions = Dict[str, int]

"""
TODO: Explain what this is...
TODO: This is actually a array in the Perl code so we could use a list
"""
SubstitutionMatrix = Dict[int, float]


def load_substitution_matrix(
    file: TextIO,
) -> Tuple[CharacterPositions, SubstitutionMatrix]:
    """
    Load the substitution matrix along with the map of character positions.

    TODO: Explain this better, like what it's used for
    TODO: Make the doctest use realistic input values

    >>> from io import StringIO
    >>> file = StringIO("a b c\\n1.0 2.0 3.0\\n4.0 5.0 6.0")
    >>> c, s = load_substitution_matrix(file)
    >>> c
    {'a': 0, 'b': 1, 'c': 2}
    >>> s
    {0: 1.0, 1: 2.0, 2: 3.0, 3: 4.0, 4: 5.0, 5: 6.0}
    """
    # First fill the character positions
    characterPositions: CharacterPositions = {}
    firstLine = next(file).replace(" ", "").replace("\t", "").strip()
    for index, character in enumerate(firstLine):
        characterPositions[character] = index

    # Then fill the substitutions
    subMatrixColumnCount = len(characterPositions)
    substitutionMatrix: SubstitutionMatrix = {}
    count = 0
    for index, line in enumerate(file):
        subScores = line.split(" ")
        for columnIndex, subScore in enumerate(subScores):
            substitutionMatrix[
                count * subMatrixColumnCount + columnIndex
            ] = float(subScore)
        count += 1

    return characterPositions, substitutionMatrix
