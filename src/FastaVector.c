#include "FastaVector.h"
#include <stdio.h>
#include <string.h>

#define FASTA_VECTOR_FILE_READ_BUFFER_SIZE 1024

// The state of the FSM used to parse the FASTA files.
enum FastaVectorFileReadState {
  FASTA_VECTOR_READ_HEADER,
  FASTA_VECTOR_READ_COMMENT,
  FASTA_VECTOR_READ_SEQUENCE,
  FASTA_VECTOR_READ_ON_NEWLINE
};

// --------------------------------
// Private function implementations
// --------------------------------

bool fastaVectorAddCharToSequenceVector(struct FastaVector *fastaVector,
                                        char c) {
  fastaVector->metadata.data[fastaVector->metadata.count - 1]
      .sequenceEndPosition++;
  return fastaVectorStringAddChar(&fastaVector->sequence, c);
}

bool fastaVectorAddCharToHeaderVector(struct FastaVector *fastaVector, char c) {
  fastaVector->metadata.data[fastaVector->metadata.count - 1]
      .headerEndPosition++;
  return fastaVectorStringAddChar(&fastaVector->header, c);
}

bool fastaVectorAddNewHeader(struct FastaVector *fastaVector) {
  return fastaVectorMetadataVectorAddMetadata(&fastaVector->metadata);
}

// -------------------------------
// Public function implementations
// -------------------------------

enum FastaVectorReturnCode fastaVectorInit(struct FastaVector *fastaVector) {
  bool result = fastaVectorStringInit(&fastaVector->sequence) &&
                fastaVectorStringInit(&fastaVector->header) &&
                fastaVectorMetadataVectorInit(&fastaVector->metadata);

  return result ? FASTA_VECTOR_OK : FASTA_VECTOR_ALLOCATION_FAIL;
}

void fastaVectorDealloc(struct FastaVector *fastaVector) {
  fastaVectorStringDealloc(&fastaVector->sequence);
  fastaVectorStringDealloc(&fastaVector->header);
  fastaVectorMetadataVectorDealloc(&fastaVector->metadata);
}

enum FastaVectorReturnCode
fastaVectorReadFasta(const char *_RESTRICT_ const fileSrc,
                     struct FastaVector *fastaVector) {

  FILE *fastaFile = fopen(fileSrc, "r");
  if (__builtin_expect(!fastaFile, false)) {
    return FASTA_VECTOR_FILE_OPEN_FAIL;
  }

  // loop over the file with a finite state machine.
  enum FastaVectorFileReadState readState = FASTA_VECTOR_READ_ON_NEWLINE;
  while (!feof(fastaFile)) {
    char charFromFile = fgetc(fastaFile);
    switch (readState) {
    case FASTA_VECTOR_READ_COMMENT:
      if (charFromFile == '\n') {
        readState = FASTA_VECTOR_READ_ON_NEWLINE;
      }
      break;

    case FASTA_VECTOR_READ_ON_NEWLINE:
      if (charFromFile == '>') {
        readState = FASTA_VECTOR_READ_HEADER;

        // add a null seperator to the sequence if this isn't the first header.
        if (fastaVector->metadata.count != 0) {
          bool addSequenceTerminatorResult =
              fastaVectorAddCharToSequenceVector(fastaVector, '\0');
          if (__builtin_expect(!addSequenceTerminatorResult, false)) {
            fclose(fastaFile);
            return FASTA_VECTOR_ALLOCATION_FAIL;
          }
        }

        bool addHeaderResult = fastaVectorAddNewHeader(fastaVector);
        if (__builtin_expect(!addHeaderResult, false)) {
          fclose(fastaFile);
          return FASTA_VECTOR_ALLOCATION_FAIL;
        }

      } else if (charFromFile == ';') {
        readState = FASTA_VECTOR_READ_COMMENT;
      } else if (__builtin_expect(charFromFile >= ' ',
                                  1)) { // only include character if they're
                                        // actually printable
        readState = FASTA_VECTOR_READ_SEQUENCE;
        bool addSequenceResult =
            fastaVectorAddCharToSequenceVector(fastaVector, charFromFile);
        if (__builtin_expect(!addSequenceResult, false)) {
          fclose(fastaFile);
          return FASTA_VECTOR_ALLOCATION_FAIL;
        }
      }
      break;

    case FASTA_VECTOR_READ_HEADER:
      if (charFromFile == ';') {
        readState = FASTA_VECTOR_READ_COMMENT;
        // since we're done with the header, null terminate it.
        bool addHeaderCharResult =
            fastaVectorAddCharToHeaderVector(fastaVector, '\0');
        if (__builtin_expect(!addHeaderCharResult, false)) {
          fclose(fastaFile);
          return FASTA_VECTOR_ALLOCATION_FAIL;
        }

      } else if (charFromFile == '\n') {
        readState = FASTA_VECTOR_READ_ON_NEWLINE;
        // since we're done with the header, null terminate it.
        bool addHeaderCharResult =
            fastaVectorAddCharToHeaderVector(fastaVector, '\0');
        if (__builtin_expect(!addHeaderCharResult, false)) {
          fclose(fastaFile);
          return FASTA_VECTOR_ALLOCATION_FAIL;
        }

      }
      // if the character looks like it belongs to the header text
      else if (__builtin_expect(charFromFile >= ' ' || charFromFile == '\t',
                                true)) {
        bool addHeaderCharResult =
            fastaVectorAddCharToHeaderVector(fastaVector, charFromFile);
        if (__builtin_expect(!addHeaderCharResult, false)) {
          fclose(fastaFile);
          return FASTA_VECTOR_ALLOCATION_FAIL;
        }
      }
      break;

    case FASTA_VECTOR_READ_SEQUENCE:
      if (charFromFile == ';') {
        readState = FASTA_VECTOR_READ_COMMENT;
      } else if (charFromFile == '\n') {
        readState = FASTA_VECTOR_READ_ON_NEWLINE;
      } else if (__builtin_expect(charFromFile >= ' ', 1)) {
        bool addSequenceCharResult =
            fastaVectorAddCharToSequenceVector(fastaVector, charFromFile);
        if (__builtin_expect(!addSequenceCharResult, false)) {
          fclose(fastaFile);
          return FASTA_VECTOR_ALLOCATION_FAIL;
        }
      }
      break;
    }
  }

  // add the final sequence terminator
  bool addSequenceTerminatorResult =
      fastaVectorAddCharToSequenceVector(fastaVector, '\0');
  if (__builtin_expect(!addSequenceTerminatorResult, false)) {
    fclose(fastaFile);
    return FASTA_VECTOR_ALLOCATION_FAIL;
  }

  fclose(fastaFile);
  return FASTA_VECTOR_OK;
}

enum FastaVectorReturnCode
fastaVectorWriteFasta(const char *_RESTRICT_ const fileSrc,
                      struct FastaVector *fastaVector,
                      uint32_t fileLineLength) {
  FILE *fastaFile = fopen(fileSrc, "w+");
  if (!fastaFile) {
    return FASTA_VECTOR_FILE_WRITE_FAIL;
  }
  for (size_t i = 0; i < fastaVector->metadata.count; i++) {
    size_t headerStartPosition =
        i == 0 ? 0 : fastaVector->metadata.data[i - 1].headerEndPosition;
    size_t headerEndPosition = fastaVector->metadata.data[i].headerEndPosition;
    size_t headerLength = headerEndPosition - headerStartPosition;

    int charWritten = fputc('>', fastaFile);
    if (charWritten != '>') {
      fclose(fastaFile);
      return FASTA_VECTOR_FILE_WRITE_FAIL;
    }

    // if the header appears to be null terminated, don't write that last
    // character (the null terminator)
    bool headerIsNullTerminated =
        fastaVector->header.charData[headerEndPosition - 1] == '\0';
    if (headerIsNullTerminated) {
      headerLength--;
    }
    size_t bytesWritten =
        fwrite(fastaVector->header.charData + headerStartPosition, sizeof(char),
               headerLength, fastaFile);
    if (bytesWritten != headerLength * sizeof(char)) {
      fclose(fastaFile);
      return FASTA_VECTOR_FILE_WRITE_FAIL;
    }
    // add the newline after the header
    fputc('\n', fastaFile);
    if (__builtin_expect(ferror(fastaFile), false)) {
      fclose(fastaFile);
      return FASTA_VECTOR_FILE_WRITE_FAIL;
    }

    // write the sequence, line by line
    const size_t sequenceStartPosition =
        i == 0 ? 0 : fastaVector->metadata.data[i - 1].sequenceEndPosition;
    // end 1 early due to sequence null seperator
    const size_t sequenceEndPosition =
        fastaVector->metadata.data[i].sequenceEndPosition - 1;

    size_t sequenceWritePosition = sequenceStartPosition;
    while (sequenceWritePosition < sequenceEndPosition) {
      size_t thisLineLength =
          (sequenceWritePosition + fileLineLength) < sequenceEndPosition
              ? fileLineLength
              : sequenceEndPosition - sequenceWritePosition;

      size_t bytesWritten =
          fwrite(fastaVector->sequence.charData + sequenceWritePosition,
                 sizeof(char), thisLineLength, fastaFile);
      if (bytesWritten != (thisLineLength * sizeof(char))) {
        fclose(fastaFile);
        return FASTA_VECTOR_FILE_WRITE_FAIL;
      }

      int charWritten = fputc('\n', fastaFile);
      if (charWritten != '\n') {
        fclose(fastaFile);
        return FASTA_VECTOR_FILE_WRITE_FAIL;
      }

      sequenceWritePosition += fileLineLength;
    }
  }
  fclose(fastaFile);
  return FASTA_VECTOR_OK;
}

enum FastaVectorReturnCode
fastaVectorAddSequenceToList(struct FastaVector *fastaVector, char *header,
                             size_t headerLength, char *sequence,
                             size_t sequenceLength) {

  // set up the metadata so it will hold the new header
  bool addHeaderResult = fastaVectorAddNewHeader(fastaVector);
  if (__builtin_expect(!addHeaderResult, false)) {
    return FASTA_VECTOR_ALLOCATION_FAIL;
  }

  // write the header
  for (size_t i = 0; i < headerLength; i++) {
    bool addHeaderCharResult =
        fastaVectorAddCharToHeaderVector(fastaVector, header[i]);
    if (__builtin_expect(!addHeaderCharResult, false)) {
      return FASTA_VECTOR_ALLOCATION_FAIL;
    }
  }

  // write the null terminator to end the header
  bool addHeaderCharResult =
      fastaVectorAddCharToHeaderVector(fastaVector, '\0');
  if (__builtin_expect(!addHeaderCharResult, false)) {
    return FASTA_VECTOR_ALLOCATION_FAIL;
  }

  // write the sequence
  for (size_t i = 0; i < sequenceLength; i++) {
    bool addSequenceCharResult =
        fastaVectorAddCharToSequenceVector(fastaVector, sequence[i]);
    if (__builtin_expect(!addSequenceCharResult, false)) {
      return FASTA_VECTOR_ALLOCATION_FAIL;
    }
  }

  // write the sequence null seperator
  bool addSequenceCharResult =
      fastaVectorAddCharToSequenceVector(fastaVector, '\0');
  if (__builtin_expect(!addSequenceCharResult, false)) {
    return FASTA_VECTOR_ALLOCATION_FAIL;
  }

  return FASTA_VECTOR_OK;
}

// return NULL for headerPtr and 0 for header length if the fastaVector does not
// have a header for the given headerIndex
void fastaVectorGetHeader(const struct FastaVector *const fastaVector,
                          size_t headerIndex, char **headerPtr,
                          size_t *headerLength) {
  if (__builtin_expect(headerIndex >= fastaVector->metadata.count, false)) {
    *headerLength = 0;
    *headerPtr = NULL;
  } else {
    const size_t headerStartPosition =
        headerIndex == 0
            ? 0
            : fastaVector->metadata.data[headerIndex - 1].headerEndPosition;
    const size_t headerEndPosition =
        fastaVector->metadata.data[headerIndex].headerEndPosition;

    // the -1 removes the terminator
    *headerLength = headerEndPosition - headerStartPosition - 1;
    *headerPtr = &fastaVector->header.charData[headerStartPosition];
  }
}

// return NULL for headerPtr and 0 for header length if the fastaVector does not
// have a header for the given headerIndex
void fastaVectorGetSequence(const struct FastaVector *const fastaVector,
                            size_t sequenceIndex, char **sequencePtr,
                            size_t *sequenceLength) {
  if (__builtin_expect(sequenceIndex >= fastaVector->metadata.count, false)) {
    *sequenceLength = 0;
    *sequencePtr = NULL;
  } else {
    const size_t sequenceStartPosition =
        sequenceIndex == 0
            ? 0
            : fastaVector->metadata.data[sequenceIndex - 1].sequenceEndPosition;
    const size_t sequenceEndPosition =
        fastaVector->metadata.data[sequenceIndex].sequenceEndPosition;

    //-1 is for the null seperator
    *sequenceLength = (sequenceEndPosition - sequenceStartPosition) - 1;
    *sequencePtr = &fastaVector->sequence.charData[sequenceStartPosition];
  }
}

bool fastaVectorGetLocalSequencePositionFromGlobal(
    const struct FastaVector *const fastaVector,
    const size_t globalSequencePosition,
    struct FastaVectorLocalPosition *localPosition) {
  if (__builtin_expect(localPosition == NULL, 0)) {
    return false;
  }
  // check to see if this is outside the bounds
  if (__builtin_expect(
          globalSequencePosition >
              fastaVector->metadata.data[fastaVector->metadata.count - 1]
                  .sequenceEndPosition,
          0)) {
    return false;
  }

  // perform binary search
  const int64_t numSequences = fastaVector->metadata.count;
  int64_t lowerBound = 0;
  int64_t upperBound = numSequences - 1;
  while (lowerBound <= upperBound) {
    const int64_t midpoint = (upperBound + lowerBound) / 2;
    const size_t midpointEndPosition =
        fastaVector->metadata.data[midpoint].sequenceEndPosition;

    if (globalSequencePosition < midpointEndPosition) {
      upperBound = midpoint - 1;
    } else {
      lowerBound = midpoint + 1;
    }
  }

  localPosition->sequenceIndex = upperBound + 1;
  localPosition->positionInSequence =
      localPosition->sequenceIndex == 0
          ? globalSequencePosition
          : globalSequencePosition -
                fastaVector->metadata.data[localPosition->sequenceIndex - 1]
                    .sequenceEndPosition;

  // if no sequence is found that contains the given position, return false to
  // show failure
  return true;
}
