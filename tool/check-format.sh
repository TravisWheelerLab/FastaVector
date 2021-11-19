#!/usr/bin/env sh

set -e

C_FILES=$(find src tests -type f -name '*.c')
H_FILES=$(find src tests -type f -name '*.h')

clang-format --dry-run --Werror ${C_FILES} ${H_FILES}
