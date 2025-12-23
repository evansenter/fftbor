// Tests for misc utility functions

#include <gtest/gtest.h>
#include <cstring>
#include "misc.h"

class MiscTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
};

TEST_F(MiscTest, XCallocAllocatesMemory) {
    // Test that xcalloc allocates memory correctly
    int* arr = (int*)xcalloc(10, sizeof(int));
    ASSERT_NE(arr, nullptr);

    // Memory should be zero-initialized
    for (int i = 0; i < 10; i++) {
        EXPECT_EQ(arr[i], 0);
    }

    free(arr);
}

TEST_F(MiscTest, XCallocLargeAllocation) {
    // Test allocation of a larger block
    double* arr = (double*)xcalloc(1000, sizeof(double));
    ASSERT_NE(arr, nullptr);

    // Memory should be zero-initialized
    for (int i = 0; i < 1000; i++) {
        EXPECT_EQ(arr[i], 0.0);
    }

    free(arr);
}

TEST_F(MiscTest, GetBasePairListSimple) {
    // Test simple hairpin structure: ((....))
    char structure[] = "((....))";
    int* bpList = getBasePairList(structure);

    ASSERT_NE(bpList, nullptr);
    EXPECT_GE(bpList[0], 0);  // No error

    // Position 1 pairs with 8 (1-indexed)
    EXPECT_EQ(bpList[1], 8);
    EXPECT_EQ(bpList[8], 1);

    // Position 2 pairs with 7
    EXPECT_EQ(bpList[2], 7);
    EXPECT_EQ(bpList[7], 2);

    // Unpaired positions (3-6) should have -1
    for (int i = 3; i <= 6; i++) {
        EXPECT_EQ(bpList[i], -1);
    }

    free(bpList);
}

TEST_F(MiscTest, GetBasePairListNested) {
    // Test nested structure: (((....)))
    char structure[] = "(((....)))";
    int* bpList = getBasePairList(structure);

    ASSERT_NE(bpList, nullptr);
    EXPECT_GE(bpList[0], 0);

    // Outermost pair: 1-10
    EXPECT_EQ(bpList[1], 10);
    EXPECT_EQ(bpList[10], 1);

    // Middle pair: 2-9
    EXPECT_EQ(bpList[2], 9);
    EXPECT_EQ(bpList[9], 2);

    // Inner pair: 3-8
    EXPECT_EQ(bpList[3], 8);
    EXPECT_EQ(bpList[8], 3);

    // Unpaired: 4-7
    for (int i = 4; i <= 7; i++) {
        EXPECT_EQ(bpList[i], -1);
    }

    free(bpList);
}

TEST_F(MiscTest, GetBasePairListMultiloop) {
    // Test multiloop structure: ((..)(..))
    char structure[] = "((..)(..))";;
    int* bpList = getBasePairList(structure);

    ASSERT_NE(bpList, nullptr);
    EXPECT_GE(bpList[0], 0);

    // Outermost pair: 1-10
    EXPECT_EQ(bpList[1], 10);
    EXPECT_EQ(bpList[10], 1);

    // First internal pair: 2-5
    EXPECT_EQ(bpList[2], 5);
    EXPECT_EQ(bpList[5], 2);

    // Second internal pair: 6-9
    EXPECT_EQ(bpList[6], 9);
    EXPECT_EQ(bpList[9], 6);

    // Unpaired positions
    EXPECT_EQ(bpList[3], -1);
    EXPECT_EQ(bpList[4], -1);
    EXPECT_EQ(bpList[7], -1);
    EXPECT_EQ(bpList[8], -1);

    free(bpList);
}

TEST_F(MiscTest, GetBasePairListAllUnpaired) {
    // Test all unpaired: ........
    char structure[] = "........";
    int n = strlen(structure);
    int* bpList = getBasePairList(structure);

    ASSERT_NE(bpList, nullptr);
    EXPECT_EQ(bpList[0], 0);  // Zero base pairs

    // All positions should be unpaired
    for (int i = 1; i <= n; i++) {
        EXPECT_EQ(bpList[i], -1);
    }

    free(bpList);
}

TEST_F(MiscTest, GetBasePairListUnbalanced) {
    // Test unbalanced structure: (((....)
    char structure[] = "(((....)";;
    int* bpList = getBasePairList(structure);

    ASSERT_NE(bpList, nullptr);
    // Should indicate an error (too many open parens)
    EXPECT_EQ(bpList[0], -2);

    free(bpList);
}

TEST_F(MiscTest, GetBasePairListTooManyClose) {
    // Test unbalanced structure: (....)))
    char structure[] = "(....)))";
    int* bpList = getBasePairList(structure);

    ASSERT_NE(bpList, nullptr);
    // Should indicate an error (too many close parens)
    EXPECT_EQ(bpList[0], -1);

    free(bpList);
}

TEST_F(MiscTest, Min2Function) {
    EXPECT_EQ(min2(5, 3), 3);
    EXPECT_EQ(min2(3, 5), 3);
    EXPECT_EQ(min2(4, 4), 4);
    EXPECT_EQ(min2(-1, 1), -1);
    EXPECT_EQ(min2(0, 0), 0);
}

TEST_F(MiscTest, Max2Function) {
    EXPECT_EQ(max2(5, 3), 5);
    EXPECT_EQ(max2(3, 5), 5);
    EXPECT_EQ(max2(4, 4), 4);
    EXPECT_EQ(max2(-1, 1), 1);
    EXPECT_EQ(max2(0, 0), 0);
}
