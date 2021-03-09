#include "../../src/FastaVector.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include "../testAssert.h"

char buffer[2048];

void readTest(char ** const testHeaderStrings, char ** const testSequenceStrings, const size_t numHeaders, const char *fastaSrc){
  struct FastaVector fastaVector;
  enum FastaVectorReturnCode returnCode = fastaVectorInit(&fastaVector);
  testAssertString(returnCode == FASTA_VECTOR_OK, "fastaInit was unsucessful for test 1");

  returnCode = fastaVectorReadFasta(fastaSrc, &fastaVector);
  testAssertString(returnCode == FASTA_VECTOR_OK, "fail on reading test1.fasta for test 1");


  for(size_t i = 0; i < numHeaders; i++){
    char *headerPtr;
    size_t headerLength;
    fastaVectorFastaGetHeader(&fastaVector, i, &headerPtr, &headerLength);
    sprintf(buffer, "header length of %zu did not match expected length %zu", headerLength, strlen(testHeaderStrings[i]));

    testAssertString(headerLength == strlen(testHeaderStrings[i]), buffer);
    bool stringsMatch = strncmp(testHeaderStrings[i], headerPtr, headerLength) == 0;
    sprintf(buffer, "header text did not match test expected value, expected '%s', got '%.*s'", testHeaderStrings[i], (int)headerLength, headerPtr);
    testAssertString(stringsMatch, buffer);

    char *sequencePtr;
    size_t sequenceLength;
    fastaVectorFastaGetSequence(&fastaVector, i, &sequencePtr, &sequenceLength);
    sprintf(buffer, "sequence length of %zu did not match expected length %zu", sequenceLength, strlen(testSequenceStrings[i]));
    testAssertString(sequenceLength == strlen(testSequenceStrings[i]), buffer);

    stringsMatch = strncmp(testSequenceStrings[i], sequencePtr, sequenceLength) == 0;
    sprintf(buffer, "sequence text did not match test expected value, expected '%s', got '%.*s'", testSequenceStrings[i], (int)sequenceLength, sequencePtr);
    testAssertString(stringsMatch, buffer);
  }

  fastaVectorDealloc(&fastaVector);
}


int main (int argc, char ** argv){
  char *testHeaderStrings[7];
  char *testSequenceStrings[7];

  printf("read test 1...\n");
  testHeaderStrings[0] = "test 1 header";
  testSequenceStrings[0] = "test1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbas";
  readTest(testHeaderStrings, testSequenceStrings, 1, "test1.fasta");

  printf("read test 2...\n");
  testHeaderStrings[0] = "test 2 header";
  testSequenceStrings[0] = "test1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbastest1sequencedataagasfnlawebrfilqawhbefrilahwbseflikhabsdlfikhbas";
  readTest(testHeaderStrings, testSequenceStrings, 1, "test2.fasta");

  printf("read test 3...\n");
  testHeaderStrings[0] = "test 2 header ";
  testSequenceStrings[0] = "awldeifbawiefbawlieskfbalkhsdbfalkhsbdflkhasbdflhjkasbdfhjkasbdfhjkasbdfiuwbeunflikzsbef";
  readTest(testHeaderStrings, testSequenceStrings, 1, "test3.fasta");

  printf("read test 4...\n");
  testHeaderStrings[0] = "test 4 header                ";
  testSequenceStrings[0] = "zxcvnzxcweryugbvmbzxcvmb";
  readTest(testHeaderStrings, testSequenceStrings, 1, "test4.fasta");

  printf("read test 5...\n");
  testHeaderStrings[0] = "header 1 ";
  testSequenceStrings[0] = "sequence1datasequence1datasequence1data";
  testHeaderStrings[1] = "  header 2 stuff with full, longer header. I'm putting some extra chars in here   1234567890-=~!@#$%^&*()_+:'\"<>?,./[]\\qwertyuiop[]\\aSDFGHJKLZXCVBNMqwertyuiop{}|ASDFGHJKKLLL:\"ZXCVBNM<>?";
  testSequenceStrings[1] = "sequence2Datasequence2Datasequence2Datasequence2Data";
  testHeaderStrings[2] = "h";
  testSequenceStrings[2] = "thisissomesequence`1234567890-=qwertyuiop[]\\\\ASDFGHJKL'ZXCVBNM,./~!@#$%^&*()_+qwertyuiop{}|asdfghjkl:\"zxcvbnm<>?some more sequence";
  readTest(testHeaderStrings, testSequenceStrings, 3, "test5.fasta");

  printf("read test 6...\n");
  testHeaderStrings[0] = "header 1";
  testSequenceStrings[0] = "sequence1DataLine1sequence1DataLine2sequence1DataLine3sequence1DataLine4sequence1DataLine5sequence1DataLine6sequence1DataLine6sequence1DataLine6sequence1DataLine6sequence1DataLine6sequence1DataLine6sequence1DataLine6";
  readTest(testHeaderStrings, testSequenceStrings, 1, "test6.fasta");


  printf("read test 7...\n");
  testHeaderStrings[0] = "header1";
  testSequenceStrings[0] = "                       seq1Dataseq1Dataseq1Data";
  testHeaderStrings[1] = "header 2";
  testSequenceStrings[1] = "seq2Dataseq2Dataseq2Dataseq2Dataseq2Dataseq2Dataseq2Data";
  testHeaderStrings[2] = "header 3";
  testSequenceStrings[2] = "seq3Dataseq3Dataseq3Dataseq3Dataseq3Dataseq3Dataseq3Dataseq3Dataseq3Dataseq3Data";
  testHeaderStrings[3] = "                  header 4 ";
  testSequenceStrings[3] = "seq4data";
  testHeaderStrings[4] = "header 5";
  testSequenceStrings[4] = "asdflkjnasdlfkjbnweriuyioy456789i8234234swedfsdfsdfpoiuytr89ol,kmjuik";
  testHeaderStrings[5] = "header6";
  testSequenceStrings[5] = "asdinhertihne";
  testHeaderStrings[6] = "header7";
  testSequenceStrings[6] = "FASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OKFASTA_VECTOR_OK";
  readTest(testHeaderStrings, testSequenceStrings, 7, "test7.fasta");

  printf("read tests completed\n");
}
