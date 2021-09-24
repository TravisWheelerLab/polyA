![PolyA CI](https://github.com/TravisWheelerLab/polyA/workflows/PolyA%20CI/badge.svg)
![GitHub](https://img.shields.io/github/license/TravisWheelerLab/polyA)
![PyPI](https://img.shields.io/pypi/v/polyA)

# AAAAAAAAAAAAAAAA (PolyA):

> **A**utomatically **A**djudicate **A**ny **A**nd **A**ll **A**rbitrary
> **A**nnotations, **A**stutely **A**djoin **A**butting **A**lignments,
> **A**nd **A**lso **A**mputate **A**nything **A**miss

A tool for adjudicating between competing annotations of biological sequences.

## About

PolyA is a software tool that adjudicates between competing alignment-based
annotations. PolyA communicates annotation certainty, identifies large scale
rearrangements, and detects boundaries between neighboring elements.

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

# or
make check
```

### Documentation

We use Sphinx for project documentation. Run `make docs` to update the
documentation sources and build the HTML. Use `make docs-serve` to serve the
documentation locally.

### Release Process

Generating a new release is a four step process. The first step MUST be
completed first, but the rest may be done in any order.

First, update the value of `VERSION` in `polyA/_version.py`, incrementing the
various portions of the version number accordingly. We would like to follow
[Semantic Versioning](https://semver.org), at least in general. Commit your
changes to the `master` branch.

Second, [create](https://github.com/TravisWheelerLab/polyA/releases/new) a
"release" on GitHub. The tag should consist of a "v", followed immediately by
the version number you chose above. For example: `v1.2.3`.

Third, run `make build-package`, followed by `make publish-package` assuming
there aren't any build errors.

Fourth, and finally, create and push a new runner image to Docker Hub. Build the
image by running `./tool/build-runner-image.sh`, and publish it with
`./tool/push-runner-image.sh`. The version tag will be set automatically and the
`latest` tag will also be updated.

## License

BSD license. See `LICENSE`.

## Authors

[Wheeler Lab](http://wheelerlab.org) at the University of Montana.

  - Kaitlin Carey
  - Audrey Shingleton
  - George Lesica
  - Jack Roddy
  - Travis Wheeler
