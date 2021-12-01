#ifndef FASTA_VECTOR_STRING_H
#define FASTA_VECTOR_STRING_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define FASTA_VECTOR_CHAR_VECTOR_DEFAULT_CAPACITY 2048

/**
 * @brief A store for FASTA data, either headers or sequences.
 */
struct FastaVectorString {
  /**
   * @brief the array of characters held by the string.
   */
  char *charData;

  /**
   * @brief the current maximum capacity of the string.
   */
  size_t capacity;

  /**
   * @brief the number of characters actually contained in the string.
   */
  size_t count;
};

/**
 * @relates FastaVectorString
 * @brief initialize the given `FastaVectorString`.
 *
 * @param vector the `FastaVectorString` to initialize
 *
 * @return true on success, false on failure
 *
 */
bool fastaVectorStringInit(struct FastaVectorString *vector);

/**
 * @relates FastaVectorString
 * @brief deinitializes the given `FastaVectorString`.
 *
 * @param vector the `FastaVectorString` to deinitialize
 *
 */
void fastaVectorStringDealloc(struct FastaVectorString *vector);

/**
 * @relates FastaVectorString
 * @brief append a character to the `FastaVectorString`, growing the string's
 * interal buffer if necessary.
 *
 * @param vector the `FastaVectorString` to modify
 * @param c the character to append
 *
 * @return true on success, false on failure
 *
 */
bool fastaVectorStringAddChar(struct FastaVectorString *vector, const char c);

#endif
