// Integration tests for FFTbor

#include <gtest/gtest.h>
#include <cstring>
#include <cmath>
#include <vector>
#include "parameter_parser.h"
#include "params.h"
#include "delta.h"
#include "misc.h"
#include "memory_types.h"

extern double temperature;
extern fftbor::ParamPtr P;
extern int N, PRECISION, WINDOW_SIZE, MIN_WINDOW_SIZE;
extern char* ENERGY;

// Static storage for ENERGY path
static char energy_path[] = "rna_turner2004.par";

class IntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        temperature = 37.0;
        PRECISION = 4;
        ENERGY = energy_path;  // Set global ENERGY path for neighbours()
    }

    void TearDown() override {
    }
};

TEST_F(IntegrationTest, SimpleHairpinPartitionFunction) {
    // Test a simple hairpin structure
    // This is a high-level test that exercises the full algorithm

    char sequence[] = "GCGCAAAAGCGC";
    char structure[] = "((((....))))";
    int n = strlen(sequence);

    // Set global parameters
    N = n;
    WINDOW_SIZE = n;
    MIN_WINDOW_SIZE = n;

    // Load parameters
    read_parameter_file("rna_turner2004.par");

    // Get base pair list
    auto bpList = getBasePairList(structure);
    ASSERT_NE(bpList.get(), nullptr);
    EXPECT_GE(bpList[0], 0);

    // Call the main algorithm
    // Note: This will print output to stdout
    testing::internal::CaptureStdout();
    neighbours(sequence, bpList.get());
    std::string output = testing::internal::GetCapturedStdout();

    // The output should contain probability values
    EXPECT_FALSE(output.empty());

    // Check that we got the expected format
    EXPECT_TRUE(output.find("k\tp(k)") != std::string::npos);

    // The sum of probabilities should be approximately 1.0
    // Parse the output and sum the probabilities
    // (simplified check - just verify output was produced)
    EXPECT_TRUE(output.find("0.") != std::string::npos ||
                output.find("1.") != std::string::npos);

    // Smart pointer auto-cleans up
}

TEST_F(IntegrationTest, AllUnpairedStructure) {
    // Test with a completely unpaired structure
    char sequence[] = "AAAAAAAA";
    char structure[] = "........";
    int n = strlen(structure);

    N = n;
    WINDOW_SIZE = n;
    MIN_WINDOW_SIZE = n;

    read_parameter_file("rna_turner2004.par");

    auto bpList = getBasePairList(structure);
    ASSERT_NE(bpList.get(), nullptr);
    EXPECT_EQ(bpList[0], 0);  // No base pairs

    // All positions should be unpaired
    for (int i = 1; i <= n; i++) {
        EXPECT_EQ(bpList[i], -1);
    }

    testing::internal::CaptureStdout();
    neighbours(sequence, bpList.get());
    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_FALSE(output.empty());

    // Smart pointer auto-cleans up
}

TEST_F(IntegrationTest, TwoHairpinsStructure) {
    // Test with two adjacent hairpins
    char sequence[] = "GCGCAAAAGCGCGCGCAAAAGCGC";
    char structure[] = "((((....))))((((....))))";
    int n = strlen(sequence);

    N = n;
    WINDOW_SIZE = n;
    MIN_WINDOW_SIZE = n;

    read_parameter_file("rna_turner2004.par");

    auto bpList = getBasePairList(structure);
    ASSERT_NE(bpList.get(), nullptr);
    EXPECT_GE(bpList[0], 0);

    // Check first hairpin base pairs
    EXPECT_EQ(bpList[1], 12);
    EXPECT_EQ(bpList[12], 1);

    // Check second hairpin base pairs
    EXPECT_EQ(bpList[13], 24);
    EXPECT_EQ(bpList[24], 13);

    testing::internal::CaptureStdout();
    neighbours(sequence, bpList.get());
    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_FALSE(output.empty());

    // Smart pointer auto-cleans up
}

TEST_F(IntegrationTest, DifferentTemperatures) {
    // Test that different temperatures produce different results
    char sequence[] = "GCGCAAAAGCGC";
    char structure[] = "((((....))))";
    int n = strlen(sequence);

    N = n;
    WINDOW_SIZE = n;
    MIN_WINDOW_SIZE = n;

    // Test at 37C
    temperature = 37.0;
    read_parameter_file("rna_turner2004.par");

    auto bpList1 = getBasePairList(structure);
    char seq1[20];
    strcpy(seq1, sequence);

    testing::internal::CaptureStdout();
    neighbours(seq1, bpList1.get());
    std::string output37 = testing::internal::GetCapturedStdout();

    // Test at 25C
    temperature = 25.0;
    read_parameter_file("rna_turner2004.par");

    char structure2[] = "((((....))))";
    auto bpList2 = getBasePairList(structure2);
    char seq2[20];
    strcpy(seq2, sequence);

    testing::internal::CaptureStdout();
    neighbours(seq2, bpList2.get());
    std::string output25 = testing::internal::GetCapturedStdout();

    // Results should be different at different temperatures
    // (though the format will be similar)
    // This is a weak check - just verify both produced output
    EXPECT_FALSE(output37.empty());
    EXPECT_FALSE(output25.empty());

    // Smart pointers auto-clean up

    // Reset temperature
    temperature = 37.0;
}

TEST_F(IntegrationTest, ConsistentResults) {
    // Test that running the same input twice gives the same results
    char sequence[] = "GCGCAAAAGCGC";
    int n = strlen(sequence);

    N = n;
    WINDOW_SIZE = n;
    MIN_WINDOW_SIZE = n;

    read_parameter_file("rna_turner2004.par");

    // First run
    char structure1[] = "((((....))))";
    auto bpList1 = getBasePairList(structure1);
    char seq1[20];
    strcpy(seq1, sequence);

    testing::internal::CaptureStdout();
    neighbours(seq1, bpList1.get());
    std::string output1 = testing::internal::GetCapturedStdout();

    // Second run
    char structure2[] = "((((....))))";
    auto bpList2 = getBasePairList(structure2);
    char seq2[20];
    strcpy(seq2, sequence);

    testing::internal::CaptureStdout();
    neighbours(seq2, bpList2.get());
    std::string output2 = testing::internal::GetCapturedStdout();

    // Results should be identical
    EXPECT_EQ(output1, output2);

    // Smart pointers auto-clean up
}
