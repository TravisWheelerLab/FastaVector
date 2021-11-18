cmake_minimum_required(VERSION 3.21)
project(
        FastaVector
        HOMEPAGE_URL https://github.com/TravisWheelerLab/FastaVector
)

set(CMAKE_C_STANDARD 17)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY build)

set(
        PUBLIC_H_FILES

        src/FastaVector.h
)

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

set_target_properties(
        fastavector_static fastavector

        PROPERTIES
        PUBLIC_HEADER ${PUBLIC_H_FILES}
)

install(TARGETS fastavector)
