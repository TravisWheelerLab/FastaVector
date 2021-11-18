#!/usr/bin/env sh

set -e

C_FILES=$(find src test -type f -name '*.c')
H_FILES=$(find src test -type f -name '*.h')

clang-format -i ${C_FILES} ${H_FILES}
