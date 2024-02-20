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
    fastaVectorFastaGetHeader(&fastaVector, i, &charPtr, &length);
    testAssertString(length == strlen(headers[i]),
                     "header lengths did not match");
    testAssertString(strcmp(headers[i], charPtr) == 0,
                     "headers did not match exactly.");
    fastaVectorFastaGetSequence(&fastaVector, i, &charPtr, &length);
    testAssertString(length == strlen(seqs[i]),
                     "sequence lengths did not match");
    testAssertString(strcmp(seqs[i], charPtr) == 0,
                     "headers did not match exactly.");
  }

  printf("test finished\n");
}