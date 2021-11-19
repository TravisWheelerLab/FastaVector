#ifndef ASSERTS_H
#define ASSERTS_H

#include <stdbool.h>
#include <string.h>

#define BEGIN_TESTS() bool __asserts_failed = false;

#define END_TESTS()                                                            \
  if (__asserts_failed) {                                                      \
    exit(1);                                                                   \
  } else {                                                                     \
    exit(0);                                                                   \
  }

#define TEST_FAILED() __asserts_failed = true;

#define ASSERT_CMP_EQ(EXP, VAL, NAME, TYPE)                                    \
  if (!(EXP == VAL)) {                                                         \
    printf("%s... failed\n"                                                    \
           "  expected: %" #TYPE "  found:    %" #TYPE,                        \
           NAME, EXP, VAL);                                                    \
    TEST_FAILED()                                                              \
  }

#define ASSERT_INT_EQ(EXP, VAL, NAME)                                          \
  if (!(EXP == VAL)) {                                                         \
    printf("%s... failed (%d != %d)\n", NAME, VAL, EXP);                       \
    TEST_FAILED()                                                              \
  }

#define ASSERT_INT_NE(EXP, VAL, NAME)                                          \
  if (EXP == VAL) {                                                            \
    printf("%s... failed (%d == %d)\n", NAME, VAL, EXP);                       \
    TEST_FAILED()                                                              \
  }

#define ASSERT_PTR_EQ(EXP, VAL, NAME)                                          \
  if (!(EXP == VAL)) {                                                         \
    printf("%s... failed (%p != %p)\n", NAME, VAL, EXP);                       \
    TEST_FAILED()                                                              \
  }

#define ASSERT_PTR_NE(EXP, VAL, NAME)                                          \
  if (EXP == VAL) {                                                            \
    printf("%s... failed (%p == %p)\n", NAME, VAL, EXP);                       \
    TEST_FAILED()                                                              \
  }

#define ASSERT_STR_EQ(EXP, VAL, NAME)                                          \
  if (strcmp(EXP, VAL) != 0) {                                                 \
    printf("%s... failed\n%s\n!=\n%s\n", NAME, VAL, EXP);                      \
    TEST_FAILED()                                                              \
  }

#define ASSERT_STR_NE(EXP, VAL, NAME)                                          \
  if (strcmp(EXP, VAL) == 0) {                                                 \
    printf("%s... failed\n%s\n==\n%s\n", NAME, VAL, EXP);                      \
    TEST_FAILED()                                                              \
  }

#endif // ASSERTS_H