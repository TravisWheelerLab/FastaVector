#include "FastaVectorString.h"

bool fastaVectorStringResize(struct FastaVectorString *vector);

bool fastaVectorStringInit(struct FastaVectorString *vector) {
  vector->capacity = FASTA_VECTOR_CHAR_VECTOR_DEFAULT_CAPACITY;
  vector->charData = malloc(vector->capacity * sizeof(char));
  vector->count = 0;
  return vector->charData != NULL;
}

void fastaVectorStringDealloc(struct FastaVectorString *vector) {
  if (vector->charData != NULL) {
    free(vector->charData);
    vector->charData = NULL;
  }
}

bool fastaVectorStringAddChar(struct FastaVectorString *vector, const char c) {
  // we check the capacity against 1 more than the count so the vector always
  // has one additional character of storage. This is useful for adding a
  // sentinel character
  // as is the case when generating an FM-index
  if (vector->capacity == (vector->count + 1)) {
    bool allocationSuccess = fastaVectorStringResize(vector);
    if (!allocationSuccess) {
      return false;
    }
  }

  vector->charData[vector->count] = c;
  vector->count++;
  return true;
}

bool fastaVectorStringResize(struct FastaVectorString *vector) {
  size_t newCapacity = vector->capacity + (vector->capacity / 2);
  vector->charData = realloc(vector->charData, newCapacity * sizeof(char));
  vector->capacity = newCapacity;

  return vector->charData != NULL;
}
