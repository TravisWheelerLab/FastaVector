#include "../src/FastaVector.h"
#include "asserts.h"
#include <stdio.h>

int main() {
  BEGIN_TESTS()

  struct FastaVector fv;
  enum FastaVectorReturnCode rc = fastaVectorInit(&fv);

  ASSERT_INT_EQ(FASTA_VECTOR_OK, 100, "init return code")

  ASSERT_PTR_NE(NULL, fv.sequence.charData, "sequence charData init")
  ASSERT_PTR_NE(NULL, fv.header.charData, "header charData init")
  ASSERT_PTR_NE(NULL, fv.metadata.data, "metdata data init")

  ASSERT_INT_EQ(0, fv.sequence.count, "sequence count init")
  ASSERT_INT_EQ(0, fv.header.count, "header count init")
  ASSERT_INT_EQ(0, fv.metadata.count, "metadata count init")

  ASSERT_INT_EQ(FASTA_VECTOR_CHAR_VECTOR_DEFAULT_CAPACITY, fv.sequence.capacity,
                "sequence capacity init")
  ASSERT_INT_EQ(FASTA_VECTOR_CHAR_VECTOR_DEFAULT_CAPACITY, fv.header.capacity,
                "header capacity init")
  ASSERT_INT_EQ(FASTA_VECTOR_SEQUENCE_METADATA_VECTOR_DEFAULT_CAPACITY,
                fv.metadata.capacity, "metadata capacity init")

  END_TESTS()
}
