from abc import ABC, abstractmethod
import numpy as np
from typing import (
    Dict,
    Generic,
    MutableMapping,
    NamedTuple,
    Optional,
    TypeVar,
    Iterator,
)

T = TypeVar("T")


class Pos(NamedTuple):
    """
    The location of an element in a sparse matrix.
    """

    row: int
    col: int

    def __str__(self) -> str:
        return f"row {self.row}, col {self.col}"


class Dim(NamedTuple):
    """
    A table dimension, bundled together for convenience.
    """

    n_rows: Optional[int]
    n_cols: Optional[int]

    def __str__(self) -> str:
        return f"{self.n_rows} rows, {self.n_cols} cols"


class SparseTable(MutableMapping[Pos, T], Generic[T], ABC):
    """
    A data table that stores its data in a sparse manner and allows
    access by row or column.

    This is an experiment to see if we can provide a clean abstraction
    over the pseudo matrices without taking an unacceptable performance
    hit.
    """

    def __str__(self) -> str:
        return f"Table: {self.dim()}"

    @abstractmethod
    def dim(self) -> Dim:
        """
        Return the logical dimensions of the table as a tuple in the
        form of `(row, col)`. The logical dimensions define the valid
        indices and don't necessarily reflect the data (if any) in the
        table.
        """
        raise NotImplementedError()

    @abstractmethod
    def redim(self, n_rows: int, n_cols: int, force: bool = False) -> None:
        """
        Alter the logical dimensions of the table so that it has the
        given number of rows and columns. If the redimension operation
        would cause data loss, this method will raise `InvalidDimension`
        unless `force` is set to `True`.
        """
        raise NotImplementedError()

    @abstractmethod
    def n_rows(self) -> Optional[int]:
        """
        Shortcut to get the number of logical rows in the table.
        """
        raise NotImplementedError()

    @abstractmethod
    def n_cols(self) -> Optional[int]:
        """
        Shortcut to get the number of logical columns in the table.
        """
        raise NotImplementedError()

    @abstractmethod
    def n_data_rows(self) -> Optional[int]:
        """
        Return the number of rows in the table that are not empty
        (contain at least one value).
        """
        raise NotImplementedError()

    @abstractmethod
    def n_data_cols(self) -> Optional[int]:
        """
        Return the number of columns in the table that are not empty
        (contain at least one value).
        """
        raise NotImplementedError()

    @abstractmethod
    def row_data_extents(self, row: int) -> (Optional[int], Optional[int]):
        """
        Return the minimum and maximum, in that order, column indices
        that contain data in the given row (inclusive).
        """
        raise NotImplementedError()

    @abstractmethod
    def col_data_extents(self, col: int) -> (Optional[int], Optional[int]):
        """
        Return the minimum and maximum, in that order, row indices
        that contain data in the given column (inclusive).
        """
        raise NotImplementedError()

    @abstractmethod
    def iter_row_indices(self, row: int) -> Iterator[int]:
        """
        Return an iterator that will yield each of the column indices
        that hold a value in the given row. No particular order is
        guaranteed.
        """
        raise NotImplementedError()

    @abstractmethod
    def iter_col_indices(self, col: int) -> Iterator[int]:
        """
        Return an iterator that will yield each of the row indices
        that hold a value in the given column. No particular order
        is guaranteed.
        """
        raise NotImplementedError()

    @abstractmethod
    def iter_row_values(self, row: int) -> Iterator[T]:
        """
        Return an iterator that will yield each of the values in the
        given row.
        """
        raise NotImplementedError()

    @abstractmethod
    def iter_col_values(self, col: int) -> Iterator[T]:
        """
        Return an iterator that will yield each of the values in the
        given column.
        """
        raise NotImplementedError()

    @abstractmethod
    def row_empty(self, row: int) -> bool:
        """
        Return `True` if the row at the given index is empty (has no
        values stored in it) and `False` otherwise.
        """
        raise NotImplementedError()

    @abstractmethod
    def col_empty(self, col: int) -> bool:
        """
        Return `True` if the column at the given index is empty (has
        no values stored in it) and `False` otherwise.
        """
        raise NotImplementedError()


class DictSparseTable(SparseTable[T]):
    """
    >>> t = DictSparseTable[float]()
    """

    _dim: Dim

    _matrix: Dict[Pos, T]

    _row_major: Dict[int, Dict[int, T]]

    _col_major: Dict[int, Dict[int, T]]

    def __init__(
        self, n_rows: Optional[int] = None, n_cols: Optional[int] = None,
    ):
        self._dim = Dim(n_rows, n_cols)
        self._matrix = {}
        self._row_major = {}
        self._col_major = {}

    def __setitem__(self, k: Pos, v: T) -> None:
        """
        >>> t = DictSparseTable[float]()
        >>> t[Pos(1, 2)] = 1.0
        >>> t._matrix[Pos(1, 2)]
        1.0
        >>> t._row_major[1][2]
        1.0
        >>> t._col_major[2][1]
        1.0
        """
        self._matrix[k] = v

        self._row_major.setdefault(k.row, {})
        self._row_major[k.row][k.col] = v

        self._col_major.setdefault(k.col, {})
        self._col_major[k.col][k.row] = v

    def __delitem__(self, k: Pos) -> None:
        """
        Remove a value from the table (implements the `del` statement).

        >>> t = DictSparseTable[float]()
        >>> t._matrix[Pos(1, 2)] = 1.0
        >>> t._row_major[1] = {2: 1.0}
        >>> t._col_major[2] = {1: 1.0}
        >>> Pos(1, 2) in t._matrix
        True
        >>> del t[Pos(1, 2)]
        >>> Pos(1, 2) in t._matrix
        False
        >>> 1 in t._row_major
        False
        >>> 2 in t._col_major
        False
        """
        if k not in self:
            return

        del self._matrix[k]

        del self._row_major[k.row][k.col]
        if len(self._row_major[k.row]) == 0:
            del self._row_major[k.row]

        del self._col_major[k.col][k.row]
        if len(self._col_major[k.col]) == 0:
            del self._col_major[k.col]

    def __getitem__(self, k: Pos) -> T:
        """
        Retrieve the value from the given position within the table
        (implements the `[]` operator).

        >>> t = DictSparseTable[float]()
        >>> t._matrix[Pos(1, 1)] = 1.0
        >>> Pos(1, 1) in t._matrix
        True
        >>> t[Pos(1, 1)]
        1.0
        """
        return self._matrix[k]

    def __len__(self) -> int:
        """
        Fetch the total number of values in the table (implements
        the `len` built-in function).

        >>> t = DictSparseTable[float]()
        >>> t._matrix[Pos(1, 1)] = 1.0
        >>> t._matrix[Pos(2, 2)] = 2.0
        >>> len(t)
        2
        """
        return len(self._matrix)

    def __iter__(self) -> Iterator[Pos]:
        """
        Return an iterator over the filled positions within the table
        (implements the `in` statement).

        >>> t = DictSparseTable[float]()
        >>> t._matrix[Pos(1, 1)] = 1.0
        >>> t._matrix[Pos(2, 2)] = 2.0
        >>> [t._matrix[p] for p in t]
        [1.0, 2.0]
        """
        return iter(self._matrix)

    def dim(self) -> Dim:
        """
        Return the logical dimensions of the table as a tuple in the
        form of `(row, col)`. The logical dimensions define the valid
        indices and don't necessarily reflect the data (if any) in the
        table.

        >>> t = DictSparseTable[float]()
        >>> t.dim()
        Dim(n_rows=None, n_cols=None)
        >>> t = DictSparseTable(1, 2)
        >>> t.dim()
        Dim(n_rows=1, n_cols=2)
        """
        return self._dim

    def redim(self, n_rows: int, n_cols: int, force: bool = False) -> None:
        """
        Alter the logical dimensions of the table so that it has the
        given number of rows and columns. If the redimension operation
        would cause data loss, this method will raise `InvalidDimension`
        unless `force` is set to `True`.

        >>> t = DictSparseTable[float]()
        >>> t._dim
        Dim(n_rows=None, n_cols=None)
        >>> t.redim(1, 2)
        >>> t._dim
        Dim(n_rows=1, n_cols=2)
        """
        # TODO: Check to see if we're going to lose data unless force=True

        self._dim = Dim(n_rows, n_cols)

    def n_rows(self) -> Optional[int]:
        """
        Shortcut to get the number of logical rows in the table.

        >>> t = DictSparseTable[float]()
        >>> t.n_rows()
        >>> t = DictSparseTable[float](1, 2)
        >>> t.n_rows()
        1
        """
        return self._dim.n_rows

    def n_cols(self) -> Optional[int]:
        """
        Shortcut to get the number of logical columns in the table.

        >>> t = DictSparseTable[float]()
        >>> t.n_cols()
        >>> t = DictSparseTable[float](1, 2)
        >>> t.n_cols()
        2
        """
        return self._dim.n_cols

    def n_data_rows(self) -> Optional[int]:
        """
        Return the number of rows in the table that are not empty
        (contain at least one value).

        >>> t = DictSparseTable[float]()
        >>> t._row_major[0] = {0: 1.0}
        >>> t.n_data_rows()
        1
        """
        return len(self._row_major)

    def n_data_cols(self) -> Optional[int]:
        """
        Return the number of columns in the table that are not empty
        (contain at least one value).

        >>> t = DictSparseTable[float]()
        >>> t._col_major[0] = {0: 1.0}
        >>> t.n_data_cols()
        1
        """
        return len(self._col_major)

    def row_data_extents(self, row: int) -> (Optional[int], Optional[int]):
        """
        Return the minimum and maximum, in that order, column indices
        that contain data in the given row (inclusive).

        >>> t = DictSparseTable[float]()
        >>> t.row_data_extents(1)
        (None, None)
        >>> t._row_major[1] = {1: 1.0, 3: 3.0}
        >>> t.row_data_extents(1)
        (1, 3)
        """
        if row not in self._row_major:
            return None, None
        keys = self._row_major[row].keys()
        return min(keys), max(keys)

    def col_data_extents(self, col: int) -> (Optional[int], Optional[int]):
        """
        Return the minimum and maximum, in that order, row indices
        that contain data in the given column (inclusive).

        >>> t = DictSparseTable[float]()
        >>> t.col_data_extents(1)
        (None, None)
        >>> t._col_major[1] = {1: 1.0, 3: 3.0}
        >>> t.col_data_extents(1)
        (1, 3)
        """
        if col not in self._col_major:
            return None, None
        keys = self._col_major[col].keys()
        return min(keys), max(keys)

    def iter_row_indices(self, row: int) -> Iterator[int]:
        """
        Return an iterator that will yield each of the column indices
        that hold a value in the given row. No particular order is
        guaranteed.

        >>> t = DictSparseTable[float]()
        >>> t._row_major[0] = {0: 1.0, 2: 3.0}
        >>> sorted(t.iter_row_indices(0))
        [0, 2]
        >>> list(t.iter_row_indices(1))
        []
        """
        if row not in self._row_major:
            return iter(())

        return iter(self._row_major[row].keys())

    def iter_col_indices(self, col: int) -> Iterator[int]:
        """
        Return an iterator that will yield each of the row indices
        that hold a value in the given column. No particular order
        is guaranteed.

        >>> t = DictSparseTable[float]()
        >>> t._col_major[0] = {0: 1.0, 2: 3.0}
        >>> sorted(t.iter_col_indices(0))
        [0, 2]
        >>> list(t.iter_col_indices(1))
        []
        """
        if col not in self._col_major:
            return iter(())

        return iter(self._col_major[col].keys())

    def iter_row_values(self, row: int) -> Iterator[T]:
        """
        Return an iterator that will yield each of the values in the
        given row.

        TODO: Yield an "empty" value for empty ones on extra param?

        >>> t = DictSparseTable[float]()
        >>> t._row_major[0] = {1: 1.0, 3: 3.0}
        >>> t._row_major
        {0: {1: 1.0, 3: 3.0}}
        >>> list(t.iter_row_values(0))
        [1.0, 3.0]
        >>> list(t.iter_row_values(1))
        []
        """
        if row not in self._row_major:
            return iter(())
        return iter(self._row_major[row].values())

    def iter_col_values(self, col: int) -> Iterator[T]:
        """
        Return an iterator that will yield each of the values in the
        given column.

        TODO: Yield an "empty" value for empty ones on extra param?

        >>> t = DictSparseTable[float]()
        >>> t._col_major[0] = {1: 1.0, 3: 3.0}
        >>> t._col_major
        {0: {1: 1.0, 3: 3.0}}
        >>> list(t.iter_col_values(0))
        [1.0, 3.0]
        >>> list(t.iter_col_values(1))
        []
        """
        if col not in self._col_major:
            return iter(())
        return iter(self._col_major[col].values())

    def row_empty(self, row: int) -> bool:
        """
        Return `True` if the row at the given index is empty (has no
        values stored in it) and `False` otherwise.

        >>> t = DictSparseTable[float]()
        >>> t[Pos(1, 1)] = 1.0
        >>> t.row_empty(0)
        True
        >>> t.row_empty(1)
        False
        """
        return row not in self._row_major

    def col_empty(self, col: int) -> bool:
        """
        Return `True` if the column at the given index is empty (has
        no values stored in it) and `False` otherwise.

        >>> t = DictSparseTable[float]()
        >>> t[Pos(1, 1)] = 1.0
        >>> t.col_empty(0)
        True
        >>> t.col_empty(1)
        False
        """
        return col not in self._col_major


class NumpySparseTable(SparseTable[T]):

    _dim: Dim

    _matrix: np.ndarray

    _filled: np.ndarray

    def __init__(
        self, n_rows: Optional[int] = None, n_cols: Optional[int] = None,
    ):
        if n_rows is None:
            self._l_rows = 2
        else:
            self._l_rows = n_rows

        if n_cols is None:
            self._l_cols = 2
        else:
            self._l_cols = n_cols

        # TODO: This will be a matrix of float64, not T
        self._matrix = np.zeros((self._l_rows, self._l_cols))
        self._filled = np.zeros((self._l_rows, self._l_cols), bool)

    def __setitem__(self, k: Pos, v: T) -> None:
        """
        >>> t = NumpySparseTable()
        >>> t._matrix
        array([[0., 0.],
               [0., 0.]])
        >>> t._filled
        array([[False, False],
               [False, False]])
        >>> t[Pos(0, 0)] = 1.0
        >>> t[Pos(1, 1)] = 1.0
        >>> t._matrix
        array([[1., 0.],
               [0., 1.]])
        >>> t._filled
        array([[ True, False],
               [False,  True]])
        """
        self._matrix[k.row, k.col] = v
        self._filled[k.row, k.col] = True

    def __delitem__(self, k: Pos) -> None:
        """
        >>> t = NumpySparseTable()
        >>> t._matrix
        array([[0., 0.],
               [0., 0.]])
        >>> t._filled
        array([[False, False],
               [False, False]])
        >>> t._matrix[0, 0] = 1.0
        >>> t._filled[0, 0] = True
        >>> t._matrix
        array([[1., 0.],
               [0., 0.]])
        >>> t._filled
        array([[ True, False],
               [False, False]])
        >>> del t[Pos(0, 0)]
        >>> t._matrix
        array([[0., 0.],
               [0., 0.]])
        >>> t._filled
        array([[False, False],
               [False, False]])
        """
        self._matrix[k.row, k.col] = 0.0
        self._filled[k.row, k.col] = False

    def __getitem__(self, k: Pos) -> Optional[T]:
        if not self._filled[k.row, k.col]:
            raise KeyError(k)
        return self._matrix[k.row, k.col]

    def __len__(self) -> int:
        return self._filled.sum()

    def __iter__(self) -> Iterator[Pos]:
        return iter(self._matrix)

    def dim(self) -> Dim:
        return Dim(self._l_rows, self._l_cols)

    def redim(self, n_rows: int, n_cols: int, force: bool = False) -> None:
        self._matrix.shape
        self._matrix.reshape((n_rows, n_cols))

    def n_rows(self) -> Optional[int]:
        pass

    def n_cols(self) -> Optional[int]:
        pass

    def n_data_rows(self) -> Optional[int]:
        pass

    def n_data_cols(self) -> Optional[int]:
        pass

    def row_data_extents(self, row: int) -> (Optional[int], Optional[int]):
        pass

    def col_data_extents(self, col: int) -> (Optional[int], Optional[int]):
        pass

    def iter_row_indices(self, row: int) -> Iterator[int]:
        pass

    def iter_col_indices(self, col: int) -> Iterator[int]:
        pass

    def iter_row_values(self, row: int) -> Iterator[T]:
        pass

    def iter_col_values(self, col: int) -> Iterator[T]:
        pass

    def row_empty(self, row: int) -> bool:
        pass

    def col_empty(self, col: int) -> bool:
        pass
