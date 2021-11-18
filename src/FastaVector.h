#ifndef FASTA_VECTOR_H
#define FASTA_VECTOR_H

/**
 * @file FastaVector.h
 * @brief The FastaVector public interface.
 *
 */

#include "FastaVectorMetadataVector.h"
#include "FastaVectorString.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/**
 * @brief Primary struct that stores data for a given FASTA file.
 *
 */
struct FastaVector {
  /**
   * @brief Storage for sequence data.
   */
  struct FastaVectorString sequence;

  /**
   * @brief Storage for header data.
   */
  struct FastaVectorString header;

  /**
   * @brief Storage for header and sequence metadata.
   */
  struct FastaVectorMetadataVector metadata;
};

/**
 * @brief Return codes used by `FastaVector` functions.
 * 
 */
enum FastaVectorReturnCode {
  /**
   * @brief The operation succeeded
   */
  FASTA_VECTOR_OK,

  /**
   * @brief The given file could not be opened
   */
  FASTA_VECTOR_FILE_OPEN_FAIL,

  /**
   * @brief The given file could not be written
   */
  FASTA_VECTOR_FILE_WRITE_FAIL,

  /**
   * @brief The given file could not be read
   */
  FASTA_VECTOR_FILE_READ_FAIL,

  /**
   * @brief Memory allocation failed
   * 
   */
  FASTA_VECTOR_ALLOCATION_FAIL
};

/**
 * @brief A position within a single sequence contained in a `FastaVector`.
 * 
 */
struct FastaVectorLocalPosition {
  /**
   * @brief The 0-based sequence index within the `FastaVector`.
   */
  size_t sequenceIndex;

  /**
   * @brief The 0-based position within the sequence.
   */
  size_t positionInSequence;
};

// ----------
// Public API
// ----------

/**
 * @relates FastaVector
 * @brief Initialize an empty `FastaVector` struct
 *
 * @param fastaVector the vector to initialize
 * @return a status code
 *
 * The struct itself must have been allocated before this function is called.
 *
 */
enum FastaVectorReturnCode fastaVectorInit(struct FastaVector *fastaVector);

/**
 * @relates FastaVector
 * @brief Deinitialize `FastaVector` member data
 *
 * @param fastaVector the vector to deinitialize
 *
 * This function does not deallocate (free) the `FastaVector` itself, only the
 * member data.
 * 
 */
void fastaVectorDealloc(struct FastaVector *fastaVector);

/**
 * @relates FastaVector
 * @brief Load the given fasta file.
 * 
 * @param fileSrc the path of the FASTA file to load
 * @param fastaVector the vector to load the fasta into
 * @param nullTerminateHeaders whether null characters should be added to headers
 * @param nullTerminateSequences whether null characters should be added to sequences
 * 
 * @return a status code
 * 
 * Loads its full contents into the initialized fasta vector. This function
 * takes any properly formed fasta file, and loads all headers and sequences.
 *
 * If `nullTerminateHeaders` is set to true, then null characters ('\0') are
 * added after each header to easily denote the end of the header (useful for
 * printing). If `nullTerminateSequences` is set to true, null characters ('\0')
 * are added after each sequence to easily separate them (useful for printing or
 * separating the sequences in an index).
 * 
 */
enum FastaVectorReturnCode
fastaVectorReadFasta(const char *restrict const fileSrc,
                     struct FastaVector *fastaVector);

/**
 * @relates FastaVector
 * @brief Write a `FastaVector` to a file.
 * 
 * @param filePath path of the file to be written
 * @param fastaVector the FASTA data to be written
 * @param fileLineLength maximum number of characters per line of sequence
 * 
 * @return a status code
 * 
 * This function will overwrite the file at the given path.
 * 
 */
enum FastaVectorReturnCode
fastaVectorWriteFasta(const char *restrict const filePath,
                      struct FastaVector *fastaVector, uint32_t fileLineLength);

/**
 * @relates FastaVector
 * @brief Add a header and sequence to the given `FastaVector`.
 * 
 * @param fastaVector the `FastaVector` to modify
 * @param header pointer to the first character in the header
 * @param headerLength length of the header, in bytes
 * @param sequence pointer to the first character in the sequence
 * @param sequenceLength length of the sequence, in bytes
 * 
 * @return a status code
 * 
 * The header should not include the leading ">" character, it will be added
 * automatically. The header and sequence will be automatically null-terminated
 * upon insertion into the `FastaVector`.
 *
 * The return value is one of:
 * 
 *   * `FASTA_VECTOR_OK`
 *   * `FASTA_VECTOR_ALLOCATION_FAIL`
 * 
 * TODO: What if the original FV was not null-terminated?
 * 
 */
enum FastaVectorReturnCode
fastaVectorAddSequenceToList(struct FastaVector *fastaVector, char *header,
                             size_t headerLength, char *sequence,
                             size_t sequenceLength);

/**
 * @relates FastaVector
 * @brief Gets a header from the given `FastaVector`.
 * 
 * @param fastaVector the `FastaVector` to read
 * @param headerIndex the index of the header to be retrieved
 * @param headerPtr pointer to store the header
 * @param headerLength pointer to store the header length
 * 
 * If the FastaVector was built with null-terminated headers, the null
 * terminator will be included in the resulting `headerPtr` and `headerLength`.
 * 
 * The index is 0-based so, for example, index 4 is the 5th header.
 * 
 * TODO: Should we return a statuc code in case the index was out of bounds?
 * 
 */
void fastaVectorFastaGetHeader(struct FastaVector *fastaVector,
                               size_t headerIndex, char **headerPtr,
                               size_t *headerLength);

/**
 * @relates FastaVector
 * @brief Gets a sequence from the given `FastaVector`.
 * 
 * @param fastaVector the `FastaVector` to read
 * @param sequenceIndex the index of the sequence to be retrieved
 * @param sequencePtr pointer to store the sequence
 * @param sequenceLength pointer to store the sequence length
 * 
 * If the FastaVector was built with null-terminated sequences, the null
 * terminator will be included in the resulting `sequencePtr` and
 * `sequenceLength`.
 * 
 * The index is 0-based so, for example, index 4 is the 5th sequence.
 * 
 */
void fastaVectorFastaGetSequence(struct FastaVector *fastaVector,
                                 size_t sequenceIndex, char **sequencePtr,
                                 size_t *sequenceLength);

/**
 * @relates FastaVector
 * @brief Gets a local sequence index and position from a `FastaVector`.
 * 
 * @param fastaVector the `FastaVector` to inspect
 * @param globalSequencePosition position across the entire sequence collection
 * @param localPosition a pointer to the struct to load the local position data into
 * 
 * @return true on success, false on failure
 * 
 * The return value will be false if the given `globalSequencePosition` is
 * outside the range of positions in the vector of sequences. In other words, if
 * `globalSequencePosition` is greater than the sum of the lengths of all the
 * sequences in the `FastaVector`, minus 1.
 * 
 */
bool fastaVectorGetLocalSequencePositionFromGlobal(
    const struct FastaVector *const fastaVector,
    const size_t globalSequencePosition,
    struct FastaVectorLocalPosition *localPosition);

/*Advanced Public Function Prototypes
 *NOTE: these functions should not be used for normal operation, but they are
 *included here for advanced users who may need to build the data structure in a
 *character-by-character basis. These functions add elements to the vectors in a
 *controlled way that grows the internal arrays if necessary. Please review the
 *implementation for fastaVectorAddSequenceToList to see an example of how to
 *use them.
 */
bool fastaVectorAddCharToSequenceVector(struct FastaVector *fastaVector,
                                        char c);
bool fastaVectorAddCharToHeaderVector(struct FastaVector *fastaVector, char c);
bool fastaVectorAddNewHeader(struct FastaVector *fastaVector);

#endif
