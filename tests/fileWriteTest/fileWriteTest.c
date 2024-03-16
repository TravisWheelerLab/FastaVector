#include "../../src/FastaVector.h"
#include "../testAssert.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

char *buffer;

void fastaVectorFileWriteTest(size_t numSequences, const char *fileSrc);
void generateRandomFastaString(char *buffer, size_t length);
void fastaVectorWriteStaticTest();

int main() {
  char fileSrcBuffer[128];
  buffer = malloc(262144 * sizeof(char));
  fastaVectorWriteStaticTest();

  for (size_t testNum = 0; testNum < 10; testNum++) {
    sprintf(fileSrcBuffer, "test%zu.fasta", testNum);
    size_t numSequences = (rand() % 10) + 1;
    fastaVectorFileWriteTest(numSequences, fileSrcBuffer);
  }

  free(buffer);
  printf("write tests finished\n");
}

void fastaVectorWriteStaticTest() {
  const char *staticFileSrc = "statictest.fasta";
  struct FastaVector fastaVector;
  enum FastaVectorReturnCode returnCode = fastaVectorInit(&fastaVector);
  testAssertString(returnCode == FASTA_VECTOR_OK,
                   "write test init did not return FASTA_VECTOR_OK");
  char *header1 = "this is header 1!";
  char *sequence1 =
      "sequence 1 is here skdvf "
      "nsldkjfnaslkdjfnsadfasdflkjabnjsedfkawbenrfkawhsbjef asjfb";
  returnCode = fastaVectorAddSequenceToList(
      &fastaVector, header1, strlen(header1), sequence1, strlen(sequence1));

  testAssertString(
      returnCode == FASTA_VECTOR_OK,
      "adding sequence to fasta vector did not return fasta vector ok.");

  char *header2 = "2";
  char *sequence2 = "wtreyeuioytrewdrfghjgfdsafghfdgbhdgv1234567890-=!@#$%^&*()"
                    "_+`qwertyuiop[]QWERTYUIOP{}asdfghjklASDFGHJKL,./<>?`";
  returnCode = fastaVectorAddSequenceToList(
      &fastaVector, header2, strlen(header2), sequence2, strlen(sequence2));
  testAssertString(returnCode == FASTA_VECTOR_OK,
                   "adding seq 2 did not return FASTA_VECTOR_OK");

  returnCode = fastaVectorWriteFasta(staticFileSrc, &fastaVector, 62);
  testAssertString(
      returnCode == FASTA_VECTOR_OK,
      "writing static fasta to file did not return FASTA_VECTOR_OK");

  struct FastaVector fastaReadVector;
  returnCode = fastaVectorInit(&fastaReadVector);
  testAssertString(returnCode == FASTA_VECTOR_OK,
                   "fasta read vector init did not return FASTA_VECTOR_OK");

  returnCode = fastaVectorReadFasta(staticFileSrc, &fastaReadVector);
  testAssertString(
      returnCode == FASTA_VECTOR_OK,
      "fasta read vector file read did not return FASTA_VECTOR_OK");

  fastaVectorDealloc(&fastaVector);
  fastaVectorDealloc(&fastaReadVector);
}

void fastaVectorFileWriteTest(const size_t numSequences, const char *fileSrc) {
  struct FastaVector fastaVector;
  enum FastaVectorReturnCode returnCode = fastaVectorInit(&fastaVector);
  testAssertString(returnCode == FASTA_VECTOR_OK,
                   "write test init did not return FASTA_VECTOR_OK");

  char **headers = malloc(numSequences * sizeof(char *));
  char **sequences = malloc(numSequences * sizeof(char *));
  uint64_t *sequenceLengths = malloc(numSequences * sizeof(uint64_t));

  // populate the fastaVector
  for (size_t sequenceNum = 0; sequenceNum < numSequences; sequenceNum++) {
    const size_t headerLength = rand() % 100 + 1;
    const size_t sequenceLength = (rand() % 100) + 1;
    headers[sequenceNum] = malloc((headerLength + 1) * sizeof(char));
    sequences[sequenceNum] = malloc((sequenceLength) * sizeof(char));
    sequenceLengths[sequenceNum] = sequenceLength;
    generateRandomFastaString(headers[sequenceNum], headerLength);
    generateRandomFastaString(sequences[sequenceNum], sequenceLength);

    // null terminate the header
    headers[sequenceNum][headerLength] = 0;

    returnCode = fastaVectorAddSequenceToList(
        &fastaVector, headers[sequenceNum], headerLength,
        sequences[sequenceNum], sequenceLength);
  }

  // check that the fastaVector returns the correct strings
  for (size_t sequenceNum = 0; sequenceNum < numSequences; sequenceNum++) {

    const char *headerPtr = headers[sequenceNum];
    const char *sequencePtr = sequences[sequenceNum];
    const size_t headerLength = strlen(headers[sequenceNum]);
    const size_t sequenceLength = sequenceLengths[sequenceNum];

    char *headerFromVector;
    char *sequenceFromVector;
    size_t headerLengthFromVector;
    size_t sequenceLengthFromVector;
    fastaVectorGetHeader(&fastaVector, sequenceNum, &headerFromVector,
                         &headerLengthFromVector);
    fastaVectorGetSequence(&fastaVector, sequenceNum, &sequenceFromVector,
                           &sequenceLengthFromVector);
    sprintf(buffer, "header was supposed to be length %zu, but got %zu.",
            headerLength, headerLengthFromVector);
    testAssertString((headerLength + 1) == headerLengthFromVector, buffer);
    sprintf(buffer, "header was supposed to match %s, but got %.*s.", headerPtr,
            (int)headerLengthFromVector, headerFromVector);
    testAssertString(strncmp(headerFromVector, headerPtr, headerLength) == 0,
                     buffer);
    sprintf(buffer, "sequence was supposed to be length %zu, but got %zu.",
            sequenceLength, sequenceLengthFromVector);
    testAssertString((sequenceLength) == sequenceLengthFromVector, buffer);
    sprintf(buffer, "sequence was supposed to match %.*s, but got %.*s.",
            (int)sequenceLength, sequencePtr, (int)sequenceLengthFromVector,
            sequenceFromVector);
    testAssertString(
        strncmp(sequenceFromVector, sequencePtr, sequenceLength) == 0, buffer);
  }

  // write the fastaVector to file
  const uint32_t fileLineLength = (rand() % 80) + 8;
  returnCode = fastaVectorWriteFasta(fileSrc, &fastaVector, fileLineLength);

  struct FastaVector fastaReadVector;
  returnCode = fastaVectorInit(&fastaReadVector);
  testAssertString(returnCode == FASTA_VECTOR_OK,
                   "fasta read vector error, write test init did not return "
                   "FASTA_VECTOR_OK");

  returnCode = fastaVectorReadFasta(fileSrc, &fastaReadVector);
  testAssertString(returnCode == FASTA_VECTOR_OK,
                   "fail on reading the written fasta file.");

  // test the capacities and counts for the 3 vectors.
  sprintf(buffer,
          "read/write fastaVectors do not agree on header char count: write "
          "%zu, read %zu",
          fastaVector.header.count, fastaReadVector.header.count);
  testAssertString(fastaVector.header.count == fastaReadVector.header.count,
                   buffer);
  if (fastaVector.header.count != fastaReadVector.header.count) {
    size_t lastHeaderLen;
    char *lastHeaderPtr;
    fastaVectorGetHeader(&fastaVector, fastaVector.metadata.count - 1,
                         &lastHeaderPtr, &lastHeaderLen);
    printf("headers did not match, \n h1: %.*s\n h2: %.*s\n hr: %.*s\n",
           (int)fastaVector.header.count, fastaVector.header.charData,
           (int)fastaReadVector.header.count, fastaReadVector.header.charData,
           (int)lastHeaderLen, lastHeaderPtr);
  }
  sprintf(buffer,
          "read/write fastaVectors do not agree on header char capacity: write "
          "%zu, read %zu",
          fastaVector.header.capacity, fastaReadVector.header.capacity);
  testAssertString(
      fastaVector.header.capacity == fastaReadVector.header.capacity, buffer);

  sprintf(buffer,
          "read/write fastaVectors do not agree on sequence char count: write "
          "%zu, read %zu",
          fastaVector.sequence.count, fastaReadVector.sequence.count);
  testAssertString(fastaVector.sequence.count == fastaReadVector.sequence.count,
                   buffer);
  sprintf(buffer,
          "read/write fastaVectors do not agree on sequence char capacity: "
          "write %zu, read %zu",
          fastaVector.sequence.capacity, fastaReadVector.sequence.capacity);
  testAssertString(fastaVector.sequence.capacity ==
                       fastaReadVector.sequence.capacity,
                   buffer);

  sprintf(buffer,
          "read/write fastaVectors do not agree on metadata count: write %zu, "
          "read %zu",
          fastaVector.metadata.count, fastaReadVector.metadata.count);
  testAssertString(fastaVector.metadata.count == fastaReadVector.metadata.count,
                   buffer);
  sprintf(buffer,
          "read/write fastaVectors do not agree on metadata capacity: write "
          "%zu, read %zu",
          fastaVector.metadata.capacity, fastaReadVector.metadata.capacity);
  testAssertString(fastaVector.metadata.capacity ==
                       fastaReadVector.metadata.capacity,
                   buffer);

  // check the header character array for identity
  for (size_t headerCharIndex = 0;
       headerCharIndex < fastaReadVector.header.count; headerCharIndex++) {
    sprintf(buffer, "header char @ pos %zu did not match, write: %c, read %c.",
            headerCharIndex, fastaVector.header.charData[headerCharIndex],
            fastaReadVector.header.charData[headerCharIndex]);
    testAssertString(fastaVector.header.charData[headerCharIndex] ==
                         fastaReadVector.header.charData[headerCharIndex],
                     buffer);
  }

  // check the sequence character array for identity
  for (size_t sequenceCharIndex = 0;
       sequenceCharIndex < fastaReadVector.sequence.count;
       sequenceCharIndex++) {
    sprintf(buffer, "header char @ pos %zu did not match, write: %c, read %c.",
            sequenceCharIndex, fastaVector.sequence.charData[sequenceCharIndex],
            fastaReadVector.sequence.charData[sequenceCharIndex]);
    testAssertString(fastaVector.sequence.charData[sequenceCharIndex] ==
                         fastaReadVector.sequence.charData[sequenceCharIndex],
                     buffer);
  }

  // check the sequence character array for identity
  for (size_t metadataIndex = 0; metadataIndex < fastaReadVector.metadata.count;
       metadataIndex++) {
    sprintf(buffer,
            "header end position @ pos %zu did not match, write: %zu, read "
            "%zu.",
            metadataIndex,
            fastaVector.metadata.data[metadataIndex].headerEndPosition,
            fastaReadVector.metadata.data[metadataIndex].headerEndPosition);
    testAssertString(
        fastaVector.metadata.data[metadataIndex].headerEndPosition ==
            fastaReadVector.metadata.data[metadataIndex].headerEndPosition,
        buffer);

    sprintf(buffer,
            "sequence end position @ pos %zu did not match, write: %zu, read "
            "%zu.",
            metadataIndex,
            fastaVector.metadata.data[metadataIndex].sequenceEndPosition,
            fastaReadVector.metadata.data[metadataIndex].sequenceEndPosition);
    testAssertString(
        fastaVector.metadata.data[metadataIndex].sequenceEndPosition ==
            fastaReadVector.metadata.data[metadataIndex].sequenceEndPosition,
        buffer);
  }

  for (size_t i = 0; i < numSequences; i++) {
    free(headers[i]);
    free(sequences[i]);
  }
  free(headers);
  free(sequences);
  free(sequenceLengths);
  fastaVectorDealloc(&fastaVector);
  fastaVectorDealloc(&fastaReadVector);
}

// generates random printable strings, replacing the comment ';' character with
// tab characters, since they're allowed too.
void generateRandomFastaString(char *buffer, size_t length) {
  for (size_t i = 0; i < length; i++) {
    char c = rand() % 92 + 32;
    if (c == ';') { // since comments could throw things off unless they're
                    // intentional, remove it if we see one.
      c = 'a';
    }
    // don't allow the first character of a line to be a GT sign.
    else if (c == '>') {
      c = '0';
    }
    buffer[i] = c;
  }
}
