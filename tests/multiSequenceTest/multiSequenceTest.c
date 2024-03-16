#include "../../src/FastaVector.h"
#include "../testAssert.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

char *headers[] = {"test 1 header\0", "test2\0", "test3\0",
                   "test4\0",         "test5\0", "z\0"};
char *seqs[] = {"abcd\0", "efgh\0", "ijk\0", "l\0", "mnop\0", "qrs\0"};

int main() {
  printf("beginning init sequence list test\n");
  struct FastaVector fastaVector;
  enum FastaVectorReturnCode returnCode = fastaVectorInit(&fastaVector);
  testAssertString(returnCode == FASTA_VECTOR_OK,
                   "fastaInit was unsucessful for test 1");

  returnCode = fastaVectorReadFasta("test.fa", &fastaVector);
  testAssertString(returnCode == FASTA_VECTOR_OK,
                   "fail on reading test1.fasta for test 1");

  for (uint8_t i = 0; i < 6; i++) {
    size_t length;
    char *charPtr;
    fastaVectorGetHeader(&fastaVector, i, &charPtr, &length);
    testAssertString(length == strlen(headers[i]),
                     "header lengths did not match");
    testAssertString(strcmp(headers[i], charPtr) == 0,
                     "headers did not match exactly.");
    fastaVectorGetSequence(&fastaVector, i, &charPtr, &length);
    testAssertString(length == strlen(seqs[i]),
                     "sequence lengths did not match");
    testAssertString(strcmp(seqs[i], charPtr) == 0,
                     "headers did not match exactly.");
  }

  // test getting the local positions
  struct FastaVectorLocalPosition lp;
  bool locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 0, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 0,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 0,
      "locate result sequence Index did not recieve expected value");

  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 3, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 3,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 0,
      "locate result sequence Index did not recieve expected value");

  // skipping the null seperator
  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 5, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 0,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 1,
      "locate result sequence Index did not recieve expected value");

  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 8, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 3,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 1,
      "locate result sequence Index did not recieve expected value");

  // skipping the null seperator
  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 10, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 0,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 2,
      "locate result sequence Index did not recieve expected value");

  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 12, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 2,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 2,
      "locate result sequence Index did not recieve expected value");

  // skipping the null seperator
  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 14, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 0,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 3,
      "locate result sequence Index did not recieve expected value");

  // skipping the null seperator
  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 16, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 0,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 4,
      "locate result sequence Index did not recieve expected value");

  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 19, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 3,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 4,
      "locate result sequence Index did not recieve expected value");

  // skipping the null seperator
  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 21, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 0,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 5,
      "locate result sequence Index did not recieve expected value");

  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 23, &lp);
  testAssertString(locateResult, "locate result returned false");
  testAssertString(lp.positionInSequence == 2,
                   "locate result position did not recieve expected value");
  testAssertString(
      lp.sequenceIndex == 5,
      "locate result sequence Index did not recieve expected value");

  locateResult =
      fastaVectorGetLocalSequencePositionFromGlobal(&fastaVector, 26, &lp);
  testAssertString(
      !locateResult,
      "locate result after sequence end did not return false, but should have");

  printf("test finished\n");
}