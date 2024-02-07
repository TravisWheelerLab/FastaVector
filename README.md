# FastaVector

A simple C library for reading/writing fasta files. This library concatenates
the headers of the fasta file into a buffer, concatenates the sequences into a
separate buffer, and keeps a metadata buffer that allows for easy indexing into
the fasta data to retrieve the headers/sequences via simple lookup table.

## Building

This library uses CMake. To build the library, use:

```
cmake .
make
```

This will produce both static and dynamic versions of the library in the
`build/` directory. To install the dynamic version, use `make install`.

The necessary headers can be found in the src/ directory. The headers,
along with the shared library can be installed using:
```
make install
```

### Legacy Build

It is also possible to build using a handrolled Makefile, without CMake. To do
this, use `Makefile_legacy`:

```
make -f Makefile_legacy
```

## Tests

The library includes a suite of unit tests. They are found in the `tests/`
directory and can be run with `make test` after the build has been run
(`make`).

The test binaries end up in the `build/` directory with the rest of the
output. You can run a suite on its own by running the corresponding binary
from the correct test directory (since they use relative paths to access
test fixtures). For example:

```
cmake .
make
cd tests/fileReadTest
../../build/fileReadTest
```

## Formatting

Code is formatted automatically using `clang-format` version 12.0.0
or later. To run the formatter, using `./tool/run-format.sh`. To check whether
formatting is required, use `./tool/check-format.sh`.

Note that formatting will occur automatically on pull requests, so manually
running the formatter is unnecessary.

## Using

For more info and the library API, please consult src/FastaVector.h and the
test suite found in the `tests/` directory.

