# FastaVector
A simple C library for reading/writing fasta files.


## Prerequisites

* A GCC-compatible compiler
* Make


To build a dynamic library (.so file), simply call
```shell
$ make
```


To build a dynamic library for Mac (.dylib),
```shell
$ make ARCH=MAC
```


To build a static library,
```shell
$ make static
```

then, to install,
```shell
$ make install
```


For more info, please consult src/fastaVector.h
