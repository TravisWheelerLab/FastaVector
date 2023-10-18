#include "FastaVectorMetadataVector.h"

bool fastaVectorMetadataVectorResize(struct FastaVectorMetadataVector *vector);

bool fastaVectorMetadataVectorInit(struct FastaVectorMetadataVector *vector) {
  vector->capacity = FASTA_VECTOR_SEQUENCE_METADATA_VECTOR_DEFAULT_CAPACITY;
  vector->data = malloc(vector->capacity * sizeof(struct FastaVectorMetadata));
  vector->count = 0;

  return vector->data != NULL;
}

bool fastaVectorMetadataVectorAddMetadata(
    struct FastaVectorMetadataVector *vector) {
  if (vector->capacity == (vector->count + 1)) {
    bool allocationSuccess = fastaVectorMetadataVectorResize(vector);
    if (!allocationSuccess) {
      return false;
    }
  }

  if (__builtin_expect(vector->count == 0, false)) {
    vector->data[0].sequenceEndPosition = 0;
    vector->data[0].headerEndPosition = 0;
  } else {
    vector->data[vector->count].sequenceEndPosition =
        vector->data[vector->count - 1].sequenceEndPosition;
    vector->data[vector->count].headerEndPosition =
        vector->data[vector->count - 1].headerEndPosition;
  }

  vector->count++;
  return true;
}

bool fastaVectorMetadataVectorResize(struct FastaVectorMetadataVector *vector) {
  size_t newCapacity = vector->capacity + (vector->capacity / 2);
  vector->data =
      realloc(vector->data, newCapacity * sizeof(struct FastaVectorMetadata));
  vector->capacity = newCapacity;

  return vector->data != NULL;
}

void fastaVectorMetadataVectorDealloc(
    struct FastaVectorMetadataVector *vector) {
  if (vector->data != NULL) {
    free(vector->data);
    vector->data = NULL;
  }
}
