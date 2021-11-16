cmake_minimum_required(VERSION 3.18)
project(
        FastaVector
        HOMEPAGE_URL https://github.com/TravisWheelerLab/FastaVector
)

set(CMAKE_C_STANDARD 17)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY build)

set(
        H_FILES

        src/FastaVector.h
        src/FastaVectorMetadataVector.h
        src/FastaVectorString.h
)
set(
        C_FILES

        src/FastaVector.c
        src/FastaVectorMetadataVector.c
        src/FastaVectorString.c
)

add_library(
        fastavector_static STATIC

        ${C_FILES}
)

add_library(
        fastavector SHARED

        ${C_FILES}
)
install(TARGETS fastavector LIBRARY)

# Runs clang-format on the code
add_custom_target(
        format
        COMMAND clang-format -i ${C_FILES} ${H_FILES}
        VERBATIM
)

# Checks whether clang-format needs to be run and sets a non-zero
# exit code if so
add_custom_target(
        check-format
        COMMAND clang-format --dry-run --Werror ${C_FILES} ${H_FILES}
        VERBATIM
)
