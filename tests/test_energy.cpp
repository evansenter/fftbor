// Tests for energy calculation functions

#include <gtest/gtest.h>
#include <cstring>
#include <cmath>
#include "parameter_parser.h"
#include "params.h"
#include "delta.h"
#include "misc.h"

extern double temperature;
extern paramT* P;

class EnergyTest : public ::testing::Test {
protected:
    void SetUp() override {
        temperature = 37.0;
        read_parameter_file("rna_turner2004.par");
        P = scale_parameters();
    }

    void TearDown() override {
        if (P) {
            free(P);
            P = nullptr;
        }
    }
};

TEST_F(EnergyTest, StackEnergiesAreNegative) {
    // Watson-Crick stacks should be stabilizing (negative)
    // CG-CG stack (type 1-1)
    EXPECT_LT(P->stack[1][1], 0);

    // GC-GC stack (type 2-2)
    EXPECT_LT(P->stack[2][2], 0);

    // AU-AU stack (type 5-5)
    EXPECT_LT(P->stack[5][5], 0);

    // UA-UA stack (type 6-6)
    EXPECT_LT(P->stack[6][6], 0);
}

TEST_F(EnergyTest, GCStacksMoreStableThanAU) {
    // GC stacks should be more stable (more negative) than AU stacks
    // CG-CG vs AU-AU
    EXPECT_LT(P->stack[1][1], P->stack[5][5]);

    // GC-GC vs AU-AU
    EXPECT_LT(P->stack[2][2], P->stack[5][5]);
}

TEST_F(EnergyTest, HairpinEnergySizeDependent) {
    // Hairpin energy should generally increase with size up to 30
    // (though there can be some variation due to special terms)

    // Very small hairpins (size 3) should be less stable than medium ones
    // Actually, size 3 has a penalty, larger is often better up to a point
    EXPECT_GT(P->hairpin[3], 0);  // All hairpins destabilizing

    // Size 4 (tetraloop) can be favorable due to special bonuses
    // but the base energy is still positive
    EXPECT_GE(P->hairpin[4], 0);

    // Larger hairpins cost more
    EXPECT_LE(P->hairpin[4], P->hairpin[10]);
    EXPECT_LE(P->hairpin[10], P->hairpin[20]);
}

TEST_F(EnergyTest, BulgeEnergySizeDependent) {
    // Bulge energy should increase with size
    EXPECT_GT(P->bulge[1], 0);  // Bulge of 1 is destabilizing
    EXPECT_LE(P->bulge[1], P->bulge[2]);
    EXPECT_LE(P->bulge[2], P->bulge[5]);
    EXPECT_LE(P->bulge[5], P->bulge[10]);
}

TEST_F(EnergyTest, InternalLoopEnergySizeDependent) {
    // Internal loop energy should increase with size
    EXPECT_GT(P->internal_loop[4], 0);  // 2x2 internal loop
    EXPECT_LE(P->internal_loop[4], P->internal_loop[8]);
    EXPECT_LE(P->internal_loop[8], P->internal_loop[16]);
}

TEST_F(EnergyTest, TerminalAUPenalty) {
    // Terminal AU penalty should be positive (destabilizing)
    EXPECT_GT(P->TerminalAU, 0);

    // It should be around 50 cal/mol (Turner parameters)
    EXPECT_GE(P->TerminalAU, 40);
    EXPECT_LE(P->TerminalAU, 60);
}

TEST_F(EnergyTest, MultiloopParameters) {
    // Multiloop closing should be positive (destabilizing)
    EXPECT_GT(P->MLclosing, 0);

    // MLbase (per unpaired nt) should be small but positive
    EXPECT_GE(P->MLbase, 0);

    // MLintern should be defined for all pair types
    for (int i = 1; i <= 6; i++) {
        // MLintern can be 0 or positive
        EXPECT_GE(P->MLintern[i], 0);
    }
}

TEST_F(EnergyTest, NinioParameters) {
    // Ninio parameters for asymmetric internal loops
    // ninio[2] is the main asymmetry parameter
    EXPECT_GT(P->ninio[2], 0);
}

TEST_F(EnergyTest, Int11Energies) {
    // 1x1 internal loops should have defined energies
    // CG-CG with A-A mismatch (indices 1,1,1,1)
    // These can be positive or negative
    int energy = P->int11[1][1][1][1];
    // Just check it's not INF or some error value
    EXPECT_LT(std::abs(energy), 10000);
}

TEST_F(EnergyTest, Int21Energies) {
    // 2x1 internal loops should have defined energies
    int energy = P->int21[1][1][1][1][1];
    EXPECT_LT(std::abs(energy), 10000);
}

TEST_F(EnergyTest, Int22Energies) {
    // 2x2 internal loops should have defined energies
    int energy = P->int22[1][1][1][1][1][1];
    EXPECT_LT(std::abs(energy), 10000);
}

TEST_F(EnergyTest, MismatchEnergies) {
    // Mismatch energies for hairpins
    // CG pair with A,A mismatch
    int mismatchH = P->mismatchH[1][1][1];
    EXPECT_LT(std::abs(mismatchH), 1000);

    // Mismatch energies for interior loops
    int mismatchI = P->mismatchI[1][1][1];
    EXPECT_LT(std::abs(mismatchI), 1000);
}

TEST_F(EnergyTest, DangleEnergies) {
    // Dangle energies (5' and 3' dangles)
    // These can be stabilizing (negative) or destabilizing (positive)
    int dangle5 = P->dangle5[1][1];  // 5' dangle on CG pair with A
    int dangle3 = P->dangle3[1][1];  // 3' dangle on CG pair with A

    EXPECT_LT(std::abs(dangle5), 1000);
    EXPECT_LT(std::abs(dangle3), 1000);
}
