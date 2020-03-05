# AAAAAAAAAAAAAAAA (PolyA)

> **A**utomatically **A**djudicate **A**ny **A**nd **A**ll **A**rbitrary
> **A**nnotations, **A**stutely > **A**djoin **A**butting **A**lignments,
> **A**nd **A**lso **A**mputate **A**nything **A**miss.

## About

**TODO:** Kaitlin

### The Algorithm

**TODO:** Kaitlin - we can also link to whatever paper(s) are published

## Using

**TODO:** Kaitlin, George

## Development

This project uses [Pipenv](https://pipenv.readthedocs.io/en/latest/). It must be installed before
any of the Makefile targets or the commands described later in this section will work.

### Makefile

A Makefile is available, run `make` or `make help` in the project root to see the available targets.

### Dependencies

If you prefer not to use the Makefile, the following commands will allow you to manage dependencies.

```
# Get to work (from within project directory)
# This drops you into a shell inside the project virtual environment
pipenv shell

# Fetch runtime dependencies
pipenv install

# Fetch runtime and development dependencies
pipenv install --dev

# Add a runtime dependency
pipenv install <package>

# Add a development dependency
pipenv install --dev package
```

### Running Perl-to-Python Tests

```
pipenv install # in the root
pipenv shell
pip install -e . # in the root
```

Then you should be able to run `RunTests.sh` in `test_inputs` successfully.

### Unit Tests

Unit tests use [pytest](https://pytest.org/en/latest/). The linked
documentation contains examples of how to use it to both write and run tests.
The [examples](https://pytest.org/en/latest/example/index.html) are quite
extensive.

Documentation tests are also fine for simple cases.

```
# Run tests
pipenv run python -m pytest

# or
make check-fast check-slow
```

## License

MIT license. See `LICENSE`.

## Authors

[Wheeler Lab](http://wheelerlab.org) at the University of Montana.

  - Travis Wheeler
  - Kaitlin Carey
  - George Lesica

