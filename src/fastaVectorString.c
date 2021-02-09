#include "fastaVectorString.h"


bool fastaVectorStringResize(struct FastaVectorString *vector);


bool fastaVectorStringInit(struct FastaVectorString *vector){
  vector->capacity = FASTA_VECTOR_CHAR_VECTOR_DEFAULT_CAPACITY;
  vector->charData = malloc(vector->capacity * sizeof(char));
  vector->count = 0;
  return vector->charData != NULL;
}

void fastaVectorStringDealloc(struct FastaVectorString *vector){
  free(vector->charData);
}


bool fastaVectorStringAddChar(struct FastaVectorString *vector, const char c){
  if(vector->capacity == vector->count){
    bool allocationSuccess = fastaVectorStringResize(vector);
    if(!allocationSuccess){
      return false;
    }
  }

  vector->charData[vector->count] = c;
  vector->count++;
  return true;
}


bool fastaVectorStringResize(struct FastaVectorString *vector){
  size_t newCapacity = vector->capacity + (vector->capacity/2);
  vector->charData = realloc(vector->charData, newCapacity* sizeof(char));
  vector->capacity = newCapacity;

  return vector->charData != NULL;
}
