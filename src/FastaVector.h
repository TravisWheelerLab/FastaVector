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

#ifdef __cplusplus
#define _RESTRICT_ __restrict__
#else
#define _RESTRICT_ restrict
#endif

/// @brief Primary struct that stores data for a given FASTA file.
struct FastaVector {
  /// @brief Storage for the fasta sequence data.
  struct FastaVectorString sequence;

  /// @brief Storage for the fasta header data.
  struct FastaVectorString header;

  /// @brief Storage for the metadata describing the FastaVector.
  struct FastaVectorMetadataVector metadata;
};

/// @brief Return codes used by `FastaVector` functions.
enum FastaVectorReturnCode {
  /// @brief The operation completed successfully
  FASTA_VECTOR_OK,

  /// @brief A given file could not be opened
  FASTA_VECTOR_FILE_OPEN_FAIL,

  /// @brief An error was encountered when attempting to write to a file.
  FASTA_VECTOR_FILE_WRITE_FAIL,

  /// @brief An error was encountered when attempting to read from a file
  FASTA_VECTOR_FILE_READ_FAIL,

  /// @brief A dynamic memory allocation was attempted, but failed.
  FASTA_VECTOR_ALLOCATION_FAIL
};

/// @brief A position within a single sequence contained in a `FastaVector`.
struct FastaVectorLocalPosition {
  /// @brief The 0-based sequence index within the `FastaVector`.
  size_t sequenceIndex;

  /// @brief The 0-based position within the sequence.
  size_t positionInSequence;
};

// ----------
// Public API
// ----------

/**
 * @relates FastaVector
 * @brief Initialize an empty `FastaVector` struct
 *
 * @param fastaVector Pointer to the vector to initialize
 * @return Possible returns are FASTA_VECTOR_OK and
 * FASTA_VECTOR_ALLOCATION_FAIL.
 *
 * The struct itself must have been allocated before this function is called.
 *
 */
enum FastaVectorReturnCode fastaVectorInit(struct FastaVector *fastaVector);

/**
 * @relates FastaVector
 * @brief Deallocates the data contained in the `FastaVector` member data.
 *
 * @param fastaVector Pointer to the FastaVector struct whose member data will
 * be deallocated.
 *
 * This function does not deallocate (free) the `FastaVector` itself, only the
 * member data. This means that if the the FastaVector struct its self was
 * dynamically allocated, the struct will still need to be freed manually.
 *
 */
void fastaVectorDealloc(struct FastaVector *fastaVector);

/**
 * @relates FastaVector
 * @brief Attempts to load the data inside the given fasta file into the given
 * FastaVector struct.
 *
 * @param fileSrc the path of the FASTA file to load
 * @param fastaVector Struct to load the FASTA data into.
 *
 * @return Possible returns are FASTA_VECTOR_OK, FASTA_VECTOR_ALLOCATION_FAIL,
 * and FASTA_VECTOR_FILE_OPEN_FAIL.
 *
 * Loads its full contents into the initialized fasta vector. This function
 * takes any properly formed fasta file, and loads all headers and sequences.
 *
 */
enum FastaVectorReturnCode
fastaVectorReadFasta(const char *_RESTRICT_ const fileSrc,
                     struct FastaVector *fastaVector);

/**
 * @relates FastaVector
 * @brief Write the data contained in a  `FastaVector` struct to the given file
 * location.
 *
 * @param filePath Path of the file to be written. This string should be
 * null-terminated.
 * @param fastaVector the FastaVector struct to be written.
 * @param fileLineLength maximum number of characters per line of sequence
 *
 * @return Possible returns are FASTA_VECTOR_OK and
 * FASTA_VECTOR_FILE_WRITE_FAIL.
 *
 * This function will overwrite the file at the given path.
 *
 */
enum FastaVectorReturnCode
fastaVectorWriteFasta(const char *_RESTRICT_ const filePath,
                      struct FastaVector *fastaVector, uint32_t fileLineLength);

/**
 * @relates FastaVector
 * @brief Add a header and sequence to the given `FastaVector`.
 *
 * @param fastaVector the `FastaVector` to append the header and sequence to.
 * @param header pointer to the first character in the header.
 * @param headerLength length of the header, in characters.
 * @param sequence pointer to the first character in the sequence.
 * @param sequenceLength length of the sequence, in characters.
 *
 * @return Possible returns are FASTA_VECTOR_OK or FASTA_VECTOR_ALLOCATION_FAIL.
 *
 * The header should not include the leading ">" character, it will be added
 * automatically. The header and sequence will be automatically null-terminated
 * upon insertion into the `FastaVector`.
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
 * @param fastaVector the `FastaVector` to read a header from.
 * @param headerIndex the index of the header to be retrieved (zero indexed).
 * @param headerPtr pointer to the char pointer that will be set to the start of
 * the header.
 * @param headerLength pointer to store the header length into.
 *
 *
 */
void fastaVectorGetHeader(const struct FastaVector *const fastaVector,
                          size_t headerIndex, char **headerPtr,
                          size_t *headerLength);

/**
 * @relates FastaVector
 * @brief Gets a sequence from the given `FastaVector`.
 *
 * @param fastaVector the `FastaVector` to read a sequence from.
 * @param sequenceIndex the index of the sequence to be retrieved
 * (zero-indexed).
 * @param sequencePtr pointer to the char pointer that will be set to the start
 * of the sequence. Ptr will be null if sequenceIndex was out of range.
 * @param sequenceLength pointer to store the sequence length into, or 0 if
 * if sequenceIndex was out of range.
 *
 * TODO: Should we return a static code in case the index was out of bounds?
 *
 */
void fastaVectorGetSequence(const struct FastaVector *const fastaVector,
                            size_t sequenceIndex, char **sequencePtr,
                            size_t *sequenceLength);

/**
 * @relates FastaVector
 * @brief Given a global position across all sequences in the FastaVector,
 * returns the sequence index and position in that sequence.
 *
 * @param fastaVector the `FastaVector` to inspect
 * @param globalSequencePosition position across the entire sequence collection.
 * @param localPosition a pointer to the struct to load the local position data
 * into
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

#endif
