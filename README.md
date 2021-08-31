# FastaVector

A simple C library for reading/writing fasta files. This library concatenates
the headers of the fasta file into a buffer, concatenates the sequences into a
separate buffer, and keeps a metadata buffer that allows for easy indexing into
the fasta data to retrieve the headers/sequences via simple lookup table.

## Prerequisites

* A GCC-compatible compiler
* Make

To build a dynamic library (.so file), simply call

```shell
$ make
```

To build a static library,
```shell
$ make static
```

then, to install,
```shell
$ make install
```

For more info and the library API, please consult src/FastaVector.h

