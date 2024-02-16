TEST_NAME = fileWriteTest
TEST_SRC	= $(TEST_NAME).c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	= -std=c11  -Wall -mtune=native -g
EXE 		= $(TEST_NAME).out

.PHONY: $(TEST_NAME)
insertTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
