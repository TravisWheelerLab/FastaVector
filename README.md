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

## Tests

The library includes a suite of unit tests. They are found in the `tests/`
directory and can be run with `make test` after the build has been run
(`make`).

### Legacy Build

It is also possible to build using a handrolled Makefile, without CMake. To do
this, use `Makefile_legacy`:

```
make -f Makefile_legacy
```

## Using

For more info and the library API, please consult src/FastaVector.h and the
test suite found in the `tests/` directory.
