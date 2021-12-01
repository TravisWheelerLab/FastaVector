#ifndef FASTA_VECTOR_METADATA_VECTOR_H
#define FASTA_VECTOR_METADATA_VECTOR_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define FASTA_VECTOR_SEQUENCE_METADATA_VECTOR_DEFAULT_CAPACITY 128

struct FastaVectorMetadata {
  size_t headerEndPosition;
  size_t sequenceEndPosition;
};

struct FastaVectorMetadataVector {
  struct FastaVectorMetadata *data;
  size_t capacity;
  size_t count;
};

/* Function:  fastaVectorMetadataVectorInit
 * --------------------
 * initializes the member data for the given metadata vector.
 *
 *  Inputs:
 *    vector: fastaVectorString to initialize.
 *
 *  Returns:
 *    true on allocation success, false on allocation failure.
 */
bool fastaVectorMetadataVectorInit(struct FastaVectorMetadataVector *vector);

/* Function:  fastaVectorMetadataVectorDealloc
 * --------------------
 * deallocates the member data for the given metadata vector.
 *
 *  Inputs:
 *    vector: fastaVectorString to deallocate.
 */
void fastaVectorMetadataVectorDealloc(struct FastaVectorMetadataVector *vector);

/* Function:  fastaVectorMetadataVectorAddMetadata
 * --------------------
 * deallocates the member data for the given metadata vector.
 *  Note that this does not populate any of the other member data. This function
 *    just updates the count of metadata to allow for population, and allocates
 * additional memory for the metadata if necessary.
 *
 *  Inputs:
 *    vector: fastaVectorString to deallocate.
 */
bool fastaVectorMetadataVectorAddMetadata(
    struct FastaVectorMetadataVector *vector);

#endif
