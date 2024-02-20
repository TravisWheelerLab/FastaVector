#include "../../src/FastaVector.h"
#include "../testAssert.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

char buffer[2048];

char randPrintableChar() {
  char c = rand() % 92 + 32;
  if (c == ';') { // since comments could throw things off unless they're
                  // intentional, remove it if we see one.
    c++;
  }
  return c;
}

void testInitSequenceList() {
  printf("beginning init sequence list test\n");
  resetTestNum();
  struct FastaVector fastaVector;

  enum FastaVectorReturnCode initResult = fastaVectorInit(&fastaVector);
  testAssertString(initResult == FASTA_VECTOR_OK,
                   "fastaVector init returned false");

  testAssertString(fastaVector.sequence.charData != NULL,
                   "sequence char data was null after init");
  testAssertString(fastaVector.header.charData != NULL,
                   "header char data was null after init");
  testAssertString(fastaVector.metadata.data != NULL,
                   "metadata was null after init");

  testAssertString(fastaVector.sequence.count == 0,
                   "sequence vector count was not null after init");
  testAssertString(fastaVector.header.count == 0,
                   "header vector count was not null after init");
  testAssertString(fastaVector.metadata.count == 0,
                   "metadata vector count was not null after init");

  testAssertString(fastaVector.sequence.capacity ==
                       FASTA_VECTOR_CHAR_VECTOR_DEFAULT_CAPACITY,
                   "sequence vector capacity was not the expected default "
                   "value after init");
  testAssertString(
      fastaVector.header.capacity == FASTA_VECTOR_CHAR_VECTOR_DEFAULT_CAPACITY,
      "header vector capacity was not the expected default value after init");
  testAssertString(fastaVector.metadata.capacity ==
                       FASTA_VECTOR_SEQUENCE_METADATA_VECTOR_DEFAULT_CAPACITY,
                   "metadata vector capacity was not the expected default "
                   "value after init");
  fastaVectorDealloc(&fastaVector);
}

void testInsertSequences() {
  printf("beginning Insert Sequence Test\n");
  resetTestNum();
  srand(time(NULL));
  struct FastaVector fastaVector;
  enum FastaVectorReturnCode initResult = fastaVectorInit(&fastaVector);
  testAssertString(initResult == FASTA_VECTOR_OK,
                   "fastaVector init returned false in insert test");

  size_t numSequences = rand() % 200 + 1;
  char *sequencePtrs[200];
  char *headerPtrs[200];
  size_t sequenceLengths[200];
  size_t headerLengths[200];

  for (size_t i = 0; i < numSequences; i++) {
    size_t headerLength = rand() % 300 + 1;
    size_t sequenceLength = rand() % 1000 + 1;
    headerPtrs[i] = malloc((headerLength * sizeof(char)) + 1);
    sequencePtrs[i] = malloc((sequenceLength * sizeof(char)) + 1);
    headerPtrs[i][0] = '>';
    for (size_t headerPos = 1; headerPos < headerLength; headerPos++) {
      headerPtrs[i][headerPos] = randPrintableChar();
    }

    for (size_t sequencePos = 0; sequencePos < sequenceLength; sequencePos++) {
      sequencePtrs[i][sequencePos] = randPrintableChar();
    }

    headerLengths[i] = headerLength;
    sequenceLengths[i] = sequenceLength;
  }

  for (size_t i = 0; i < numSequences; i++) {
    enum FastaVectorReturnCode returnCode = fastaVectorAddSequenceToList(
        &fastaVector, headerPtrs[i], headerLengths[i], sequencePtrs[i],
        sequenceLengths[i]);

    sprintf(buffer,
            "adding sequence to list for seq %zu did not return "
            "FASTA_VECTOR_OK. \nheader: %.*s, \n seq: %.*s",
            i, (int)headerLengths[i], headerPtrs[i], (int)sequenceLengths[i],
            sequencePtrs[i]);
    testAssertString(returnCode == FASTA_VECTOR_OK, buffer);

    sprintf(buffer, "expected %zu items in metadata list, got %zu", i,
            fastaVector.metadata.count);
    testAssertString(fastaVector.metadata.count == (i + 1), buffer);
  }

  // check all the lengths against the metadata, and the
  for (size_t i = 0; i < numSequences; i++) {
    size_t sequenceStartPosition =
        i == 0 ? 0 : fastaVector.metadata.data[i - 1].sequenceEndPosition;
    size_t sequenceEndPosition =
        fastaVector.metadata.data[i].sequenceEndPosition - 1;
    size_t sequenceLength = sequenceEndPosition - sequenceStartPosition;

    size_t headerStartPosition =
        i == 0 ? 0 : fastaVector.metadata.data[i - 1].headerEndPosition;
    size_t headerEndPosition = fastaVector.metadata.data[i].headerEndPosition;
    size_t headerLength = headerEndPosition - headerStartPosition;

    sprintf(buffer,
            "header %zu had length %zu, but range [%zu, %zu] makes stored "
            "length %zu",
            i, headerLengths[i], headerStartPosition, headerEndPosition,
            headerLength);
    testAssertString((headerLengths[i] + 1) == headerLength, buffer);

    sprintf(buffer,
            "sequence %zu had length %zu, but range [%zu, %zu] makes stored "
            "length %zu",
            i, sequenceLengths[i], sequenceStartPosition, sequenceEndPosition,
            sequenceLength);
    testAssertString(sequenceLengths[i] == sequenceLength, buffer);

    for (size_t i = 0; i < numSequences; i++) {
      size_t sequenceStartPosition =
          i == 0 ? 0 : fastaVector.metadata.data[i - 1].sequenceEndPosition;
      char *mainSequencePtr = fastaVector.sequence.charData;
      sprintf(buffer, "sequence number %zu did not match what was stored.", i);
      bool stringsMatch =
          strncmp(sequencePtrs[i], mainSequencePtr + sequenceStartPosition,
                  sequenceLengths[i]) == 0;
      testAssertString(stringsMatch, buffer);

      size_t headerStartPosition =
          i == 0 ? 0 : fastaVector.metadata.data[i - 1].headerEndPosition;
      char *mainHeaderPtr = fastaVector.header.charData;
      sprintf(buffer, "header number %zu did not match what was stored.", i);
      stringsMatch = strncmp(headerPtrs[i], mainHeaderPtr + headerStartPosition,
                             headerLengths[i]) == 0;
      testAssertString(stringsMatch, buffer);
    }
  }

  // cleanup and dealloc
  for (size_t i = 0; i < numSequences; i++) {
    free(headerPtrs[i]);
    free(sequencePtrs[i]);
  }

  fastaVectorDealloc(&fastaVector);
}

int main() {
  printf("beginning tests\n");
  testInitSequenceList();
  testInsertSequences();
  printf("insert tests finished\n");
}
