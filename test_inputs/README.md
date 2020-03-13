# Perl to Python Validation Tests

The tests in this directory attempt to validate that the new Python version of the adjudication
code, written to mimic the Perl version, produces the same results as the Perl version of the code.
To run then, just run the `RunTests.sh` script within the PolyA virtual environment. The complete
procedure is below:

```
cd polyA
poetry install
cd test_inputs
./RunTests.sh
```

Diffs between the two versions will be saved to files ending in `.diff`. The "left" side of the diff
is the Perl version of the output, and the "right" side of the diff is Python.

