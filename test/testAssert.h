#ifndef CIMPLE_FASTA_TEST_ASSERT_H
#define CIMPLE_FASTA_TEST_ASSERT_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

size_t testNum = 0;

void testAssert(bool testSuccess) {
  if (!testSuccess) {
    printf("TEST FAIL: Test %zu returned false.\n", testNum);
  }
  testNum++;
}

void testAssertString(bool testSuccess, char *message) {
  if (!testSuccess) {
    printf("TEST FAIL: Test %zu, %s\n", testNum, message);
  }
  testNum++;
}

void resetTestNum() { testNum = 0; }

#endif
