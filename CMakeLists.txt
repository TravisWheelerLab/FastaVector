cmake_minimum_required(VERSION 3.21)
project(
        FastaVector
        HOMEPAGE_URL https://github.com/TravisWheelerLab/FastaVector
)

set(CMAKE_C_STANDARD 17)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY build)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY build)

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

# Specify different sets of header files for shared and static libraries
set_target_properties(
        fastavector
        PROPERTIES
        PUBLIC_HEADER src/FastaVector.h
)

# set_target_properties(
#         fastavector_static
#         PROPERTIES
#         PUBLIC_HEADER src/FastaVector.h
# )

install(TARGETS fastavector fastavector_static)

# Install only the necessary header files for each library
install(FILES src/FastaVectorMetadataVector.h
    DESTINATION include
)

install(FILES src/FastaVectorString.h
    DESTINATION include
)

# Add -O3 and -fPIC flags for the shared library
target_compile_options(
        fastavector
        PRIVATE
        -O3
        -fPIC
)
# Add -O3 and -fPIC flags for static library
target_compile_options(
        fastavector_static 
        PRIVATE
        -O3
        -fPIC
)



# ----------
# Unit Tests
# ----------

include(CTest)
list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")

if (BUILD_TESTING)
        # fileReadTest

        add_executable(
                fileReadTest

                tests/fileReadTest/fileReadTest.c
                ${C_FILES}
        )
        target_compile_options(
                fileReadTest

                PRIVATE
                -mtune=native
                -Wall
                -Wextra
                -g
        )
        target_link_options(
                fileReadTest

                PRIVATE
        )
        add_test(NAME fileReadTest COMMAND fileReadTest)
        set_tests_properties(
                fileReadTest

                PROPERTIES
                WORKING_DIRECTORY tests/fileReadTest
        )

        # fileWriteTest

        add_executable(
                fileWriteTest

                tests/fileWriteTest/fileWriteTest.c
                ${C_FILES}
        )
        target_compile_options(
                fileWriteTest

                PRIVATE
                -mtune=native
                -Wall
                -Wextra
                -g
        )
        target_link_options(
                fileWriteTest

                PRIVATE
        )
        add_test(NAME fileWriteTest COMMAND fileWriteTest)
        set_tests_properties(
                fileWriteTest

                PROPERTIES
                WORKING_DIRECTORY tests/fileWriteTest
        )

        # insertTest

        add_executable(
                insertTest

                tests/insertTest/insertTest.c
                ${C_FILES}
        )
        target_compile_options(
                insertTest

                PRIVATE
                -mtune=native
                -Wall
                -Wextra
                -g
        )
        target_link_options(
                insertTest

                PRIVATE
        )
        add_test(NAME insertTest COMMAND insertTest)
        set_tests_properties(
                insertTest

                PROPERTIES
                WORKING_DIRECTORY tests/insertTest
        )

        # multiSequenceTest

        add_executable(
                multiSequenceTest

                tests/multiSequenceTest/multiSequenceTest.c
                ${C_FILES}
        )
        target_compile_options(
                multiSequenceTest

                PRIVATE
                -mtune=native
                -Wall
                -Wextra
                -g
        )
        target_link_options(
                multiSequenceTest

                PRIVATE
        )
        add_test(NAME multiSequenceTest COMMAND multiSequenceTest)
        set_tests_properties(
                multiSequenceTest

                PROPERTIES
                WORKING_DIRECTORY tests/multiSequenceTest
        )
endif ()
