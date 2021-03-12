#ifndef FASTA_VECTOR_H
#define FASTA_VECTOR_H

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include "FastaVectorString.h"
#include "FastaVectorMetadataVector.h"


//main struct to store the fasta data.
struct FastaVector{
  struct FastaVectorString sequence;
  struct FastaVectorString header;
  struct FastaVectorMetadataVector metadata;
};

enum FastaVectorReturnCode{
  FASTA_VECTOR_OK, FASTA_VECTOR_FILE_OPEN_FAIL, FASTA_VECTOR_FILE_WRITE_FAIL,
  FASTA_VECTOR_FILE_READ_FAIL, FASTA_VECTOR_ALLOCATION_FAIL};

struct FastaVectorLocalPosition{
  size_t sequenceIndex;
  size_t positionInSequence;
};

/*Public API Function Prototypes*/
/* Function:  fastaVectorInit
 * --------------------
 * Initializes the given allocated (i.e., not null ptr) fastaVector.
 *
 *  Inputs:
 *    fastaVector: vector to initialize
 *
 *  Returns:
 *    FASTA_VECTOR_OK on initialization success
 *    FASTA_VECTOR_ALLOCATION_FAIL: on failure due to a failed malloc allocation.
 */
enum FastaVectorReturnCode fastaVectorInit(struct FastaVector *fastaVector);


/* Function:  fastaVectorDealloc
 * --------------------
 * deallocates the given fasta vector's member data.
 *
 *  Inputs:
 *    fastaVector: vector to deallocate
 */
void fastaVectorDealloc(struct FastaVector *fastaVector);


/* Function:  fastaVectorReadFasta
 * --------------------
 * Reads the given fasta, and loads its full contents into the initialized fasta vector.
 *  This function takes any properly formed fasta file, and loads all headers and sequences.
 *
 *  Inputs:
 *    fileSrc: string representing the location of the fasta file to load.
 *    fastaVector: vector to load the fasta into.
 *    nullTerminateHeaders: if this is set to true, null-terminating '\0's are added
 *      after each header to easily denote the end of the header (useful for printing).
 *    nullTerminateSequences: if this is set to true, null-terminating '\0's are added
 *      after each sequence to easily separate them, useful for printing or separating
 *      the sequences in an index (if terminators are replaced with ambiguity characters.)
 *
 *  Returns:
 *    FASTA_VECTOR_OK on fasta read success
 *    FASTA_VECTOR_ALLOCATION_FAIL: on failure due to a failed malloc allocation.
 *    FASTA_VECTOR_FILE_OPEN_FAIL: on failure to open the file (was the fileSrc correct?).
 */
enum FastaVectorReturnCode fastaVectorReadFasta(const char *restrict const fileSrc,
  struct FastaVector *fastaVector, const bool nullTerminateHeaders, const bool nullTerminateSequences);


/* Function:  fastaVectorWriteFasta
 * --------------------
 * takes the given, populated fastaVector, and writes it to a file.
 *  This function will overwrite the file at this location.
 *
 *  Inputs:
 *    fileSrc: string representing the location of the fasta file to write to.
 *    fastaVector: fastaVector that will be written to a file.
 *    fileLineLength: length of the lines in the file to be written. This does not apply
 *      to headers, since they cannot be multiple lines.
 *
 *  Returns:
 *    FASTA_VECTOR_OK on fasta write success
 *    FASTA_VECTOR_FILE_WRITE_FAIL: on failure to open the file for writing.
 */
enum FastaVectorReturnCode fastaVectorWriteFasta(const char *restrict const fileSrc,
  struct FastaVector *fastaVector, uint32_t fileLineLength);


/* Function:  fastaVectorAddSequenceToList
 * --------------------
 * Adds a header and sequence to the given fastaVector.
 *
 *  Inputs:
 *    fastaVector: fastaVector to add the header and sequence to.
 *    header: pointer to the first character in the header. This does not need to be null-terminated.
 *      This header string should not include the starting '>' character, it will be
 *      added automatically if written to file.
 *    headerLength: length of the header, in characters.
 *    sequence: pointer to the first character in the sequence. This does not need to be null-terminated.
 *    sequenceLength: length of the sequence, in characters.
 *
 *  Returns:
 *    FASTA_VECTOR_OK: on success.
 *    FASTA_VECTOR_ALLOCATION_FAIL: on failure to allocate the memory to grow the internal data vectors.
 */
enum FastaVectorReturnCode fastaVectorAddSequenceToList(struct FastaVector *fastaVector, char *header,
  size_t headerLength, char *sequence, size_t sequenceLength);


/* Function:  fastaVectorFastaGetHeader
 * --------------------
 * Gets a header from the given fastaVector. The headerIndex argument determines which header is retrieved.
 *  If the FastaVector was built with null-terminated headers, the null terminator will be included in the
 *  returned headerPtr and headerLength.
 *
 *  Inputs:
 *    fastaVector: fastaVector extract a header from.
 *    headerIndex: the index of the header to extract, e.g., a header index of 4 retrieves the
 *      5th header in the fastaVector.
 *    headerPtr: pointer to the char pointer that will be set to the start of the header.
 *    headerLength: pointer to the size_t variable that will be set to the header's length
 */
void fastaVectorFastaGetHeader(struct FastaVector *fastaVector, size_t headerIndex, char **headerPtr, size_t *headerLength);


/* Function:  fastaVectorFastaGetSequence
 * --------------------
 * Gets a sequence from the given fastaVector. The headerIndex argument determines which sequence is retrieved.
 *  If the FastaVector was built with null-terminated sequences, the null terminator will be included in the
 *  returned sequencePtr and sequenceLength.
 *
 *  Inputs:
 *    fastaVector: fastaVector extract a sequence from.
 *    headerIndex: the index of the sequence to extract, e.g., a sequence index of 4 retrieves the
 *      5th sequence in the fastaVector.
 *    sequencePtr: pointer to the char pointer that will be set to the start of the sequence.
 *    sequenceLength: pointer to the size_t variable that will be set to the sequence's length
 */
void fastaVectorFastaGetSequence(struct FastaVector *fastaVector, size_t sequenceIndex, char **sequencePtr, size_t *sequenceLength);


/* Function:  fastaVectorGetLocalSequencePositionFromGlobal
 * --------------------
 * Gets the sequence index, and position in that sequence, that corresponds to
 *  the given globalSequencePosition in a full fastaVector.
 *
 *  Inputs:
 *    fastaVector: fastaVector to that the position looks into.
 *    globalSequencePosition: position across the entire sequence collection.
 *    localPosition:  ptr to the struct to load the local position data into.
 *
 *  Returns:
 *    true on success, false if the given globalSequencePosition is outside the range of positions
 *      in the vector of sequences (i.e., globalSequencePosition >= the combined length of all the sequences).
 */
bool fastaVectorGetLocalSequencePositionFromGlobal(
  const struct FastaVector *const fastaVector, const size_t globalSequencePosition, struct FastaVectorLocalPosition *localPosition);

/*Advanced Public Function Prototypes
*NOTE: these functions should not be used for normal operation, but they are included here
* for advanced users who may need to build the data structure in a character-by-character basis.
* These functions add elements to the vectors in a controlled way that grows the internal arrays if necessary.
* Please review the implementation for fastaVectorAddSequenceToList to see an example of how to use them.
*/
bool fastaVectorAddCharToSequenceVector(struct FastaVector *fastaVector, char c);
bool fastaVectorAddCharToHeaderVector(struct FastaVector *fastaVector, char c);
bool fastaVectorAddNewHeader(struct FastaVector *fastaVector);



#endif
