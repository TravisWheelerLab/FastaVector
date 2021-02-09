#ifndef FASTA_VECTOR_STRING_H
#define FASTA_VECTOR_STRING_H

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>


#define FASTA_VECTOR_CHAR_VECTOR_DEFAULT_CAPACITY 2048

struct FastaVectorString{
  char *charData;
  size_t capacity;
  size_t count;

};

/* Function:  fastaVectorStringInit
 * --------------------
 * initializes the data vector for the given fastaVectorString
 *
 *  Inputs:
 *    vector: fastaVectorString to initialize.
 *    header: pointer to the first character in the header. This does not need to be null-terminated.
 *    headerLength: length of the header, in characters.
 *    sequence: pointer to the first character in the sequence. This does not need to be null-terminated.
 *    sequenceLength: length of the sequence, in characters.
 *
 *  Returns:
 *    true on allocation success, false on allocation failure.
 */
bool fastaVectorStringInit(struct FastaVectorString *vector);


/* Function:  fastaVectorStringDealloc
 * --------------------
 * deallocates the member data for the vectorString.
 *
 *  Inputs:
 *    vector: fastaVectorString to deallocate.
 */
void fastaVectorStringDealloc(struct FastaVectorString *vector);


/* Function:  fastaVectorStringAddChar
 * --------------------
 * appends a character to the fastaVectorString. This will grow the string's interal buffer if necessary.
 *
 *  Inputs:
 *    vector: fastaVectorString to append a character to.
 *    c:      character to append.
 *
 *  Returns:
 *    true on allocation success, false on allocation failure.
 */
bool fastaVectorStringAddChar(struct FastaVectorString *vector, const char c);


#endif
