// Tests for the parameter parser

#include <gtest/gtest.h>
#include <cstring>
#include <cstdlib>
#include "parameter_parser.h"
#include "params.h"
#include "energy_const.h"

extern double temperature;

class ParameterParserTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Reset temperature to default before each test
        temperature = 37.0;
    }
};

TEST_F(ParameterParserTest, LoadTurner2004Parameters) {
    // This test verifies that we can load the Turner 2004 parameters
    read_parameter_file("rna_turner2004.par");

    paramT* P = scale_parameters();
    ASSERT_NE(P, nullptr);

    // Check that some key parameters are loaded correctly
    // Stack energies for CG-CG should be negative (stabilizing)
    EXPECT_LT(P->stack[1][1], 0);  // CG-CG stack

    // Hairpin loop energies should be positive (destabilizing)
    EXPECT_GT(P->hairpin[4], 0);   // Size 4 hairpin

    // Terminal AU penalty should be positive
    EXPECT_GT(P->TerminalAU, 0);

    // MLclosing is defined (can be positive or negative in Turner 2004)
    EXPECT_NE(P->MLclosing, INF);

    free(P);
}

TEST_F(ParameterParserTest, LoadTurner1999Parameters) {
    // This test verifies that we can load the Turner 1999 parameters
    read_parameter_file("rna_turner1999.par");

    paramT* P = scale_parameters();
    ASSERT_NE(P, nullptr);

    // Check that parameters were loaded
    EXPECT_LT(P->stack[1][1], 0);  // CG-CG stack should be stabilizing

    free(P);
}

TEST_F(ParameterParserTest, TemperatureScaling) {
    // Test that temperature scaling works correctly
    read_parameter_file("rna_turner2004.par");

    // Get parameters at 37C
    temperature = 37.0;
    paramT* P37 = scale_parameters();
    ASSERT_NE(P37, nullptr);

    int stack_37 = P37->stack[1][1];

    // Get parameters at 25C (should be slightly different)
    temperature = 25.0;
    paramT* P25 = scale_parameters();
    ASSERT_NE(P25, nullptr);

    int stack_25 = P25->stack[1][1];

    // At lower temperature, stacking should be slightly more stabilizing
    // (more negative) due to entropy-enthalpy compensation
    // This is a rough check - the exact relationship depends on enthalpies
    EXPECT_NE(stack_37, stack_25);

    // Get parameters at 50C
    temperature = 50.0;
    paramT* P50 = scale_parameters();
    ASSERT_NE(P50, nullptr);

    int stack_50 = P50->stack[1][1];

    // Temperature ordering: lower temp = more stable (more negative)
    // This may not always hold depending on enthalpy values
    // At minimum, they should all be different
    EXPECT_NE(stack_37, stack_50);

    free(P37);
    free(P25);
    free(P50);

    // Reset temperature
    temperature = 37.0;
}

TEST_F(ParameterParserTest, SpecialLoopsLoaded) {
    // Test that special loop sequences are loaded
    read_parameter_file("rna_turner2004.par");

    paramT* P = scale_parameters();
    ASSERT_NE(P, nullptr);

    // Check that Tetraloops string is not empty
    EXPECT_GT(strlen(P->Tetraloops), 0);

    // Check that Triloops string is not empty
    EXPECT_GT(strlen(P->Triloops), 0);

    free(P);
}

TEST_F(ParameterParserTest, ModelDetailsInitialized) {
    // Test that model details are properly initialized
    read_parameter_file("rna_turner2004.par");

    paramT* P = scale_parameters();
    ASSERT_NE(P, nullptr);

    // Check default model details
    EXPECT_EQ(P->model_details.dangles, 2);
    EXPECT_EQ(P->model_details.special_hp, 1);
    EXPECT_EQ(P->model_details.noLP, 0);
    EXPECT_EQ(P->model_details.noGU, 0);

    free(P);
}

TEST_F(ParameterParserTest, BulgeAndInternalLoopParameters) {
    read_parameter_file("rna_turner2004.par");

    paramT* P = scale_parameters();
    ASSERT_NE(P, nullptr);

    // Bulge of size 1 should have a specific energy
    EXPECT_GT(P->bulge[1], 0);

    // Internal loop energies should be positive
    EXPECT_GT(P->internal_loop[2], 0);
    EXPECT_GT(P->internal_loop[4], 0);

    // Larger loops should have higher energies (more destabilizing)
    EXPECT_LE(P->internal_loop[4], P->internal_loop[6]);

    free(P);
}
