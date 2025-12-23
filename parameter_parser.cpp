#include "parameter_parser.h"
#include "energy_const.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cctype>

// Global storage for raw (unscaled) energy parameters at 37C
static struct {
    int stack[NBPAIRS+1][NBPAIRS+1];
    int stack_enthalpies[NBPAIRS+1][NBPAIRS+1];
    int hairpin[31];
    int hairpin_enthalpies[31];
    int bulge[MAXLOOP+1];
    int bulge_enthalpies[MAXLOOP+1];
    int internal_loop[MAXLOOP+1];
    int internal_loop_enthalpies[MAXLOOP+1];
    int mismatchH[NBPAIRS+1][5][5];
    int mismatchH_enthalpies[NBPAIRS+1][5][5];
    int mismatchI[NBPAIRS+1][5][5];
    int mismatchI_enthalpies[NBPAIRS+1][5][5];
    int mismatch1nI[NBPAIRS+1][5][5];
    int mismatch1nI_enthalpies[NBPAIRS+1][5][5];
    int mismatch23I[NBPAIRS+1][5][5];
    int mismatch23I_enthalpies[NBPAIRS+1][5][5];
    int mismatchM[NBPAIRS+1][5][5];
    int mismatchM_enthalpies[NBPAIRS+1][5][5];
    int mismatchExt[NBPAIRS+1][5][5];
    int mismatchExt_enthalpies[NBPAIRS+1][5][5];
    int dangle5[NBPAIRS+1][5];
    int dangle5_enthalpies[NBPAIRS+1][5];
    int dangle3[NBPAIRS+1][5];
    int dangle3_enthalpies[NBPAIRS+1][5];
    int int11[NBPAIRS+1][NBPAIRS+1][5][5];
    int int11_enthalpies[NBPAIRS+1][NBPAIRS+1][5][5];
    int int21[NBPAIRS+1][NBPAIRS+1][5][5][5];
    int int21_enthalpies[NBPAIRS+1][NBPAIRS+1][5][5][5];
    int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
    int int22_enthalpies[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
    int ninio[5];
    int ninio_enthalpies[5];
    double lxc;
    int MLbase;
    int MLbase_enthalpies;
    int MLintern[NBPAIRS+1];
    int MLintern_enthalpies[NBPAIRS+1];
    int MLclosing;
    int MLclosing_enthalpies;
    int TerminalAU;
    int TerminalAU_enthalpies;
    int DuplexInit;
    int DuplexInit_enthalpies;
    int Tetraloop_E[200];
    int Tetraloop_E_enthalpies[200];
    char Tetraloops[1401];
    int num_tetraloops;
    int Triloop_E[40];
    int Triloop_E_enthalpies[40];
    char Triloops[241];
    int num_triloops;
    int Hexaloop_E[40];
    int Hexaloop_E_enthalpies[40];
    char Hexaloops[1801];
    int num_hexaloops;
    int TripleC;
    int MultipleCA;
    int MultipleCB;
    int loaded;
} raw_params;

extern double temperature;

// Helper: skip whitespace and comments
static char* skip_whitespace(char* p) {
    while (*p && isspace((unsigned char)*p)) p++;
    return p;
}

// Helper: read a line, skipping comments
static int read_line(FILE* fp, char* buf, int size) {
    while (fgets(buf, size, fp)) {
        char* p = skip_whitespace(buf);
        // Skip empty lines and comment-only lines (but not section headers)
        if (*p == '\0') continue;
        if (*p == '#' && p[1] == '#') continue; // Skip ## comments
        return 1;
    }
    return 0;
}

// Helper: parse integers from a line (handles INF values)
static int parse_int_line(char* line, int* values, int max_values) {
    int count = 0;
    char* p = line;
    while (*p && count < max_values) {
        while (*p && (isspace((unsigned char)*p) || *p == ',')) p++;
        if (*p == '\0' || *p == '/' || *p == '#') break;

        // Check for INF (case insensitive)
        if ((p[0] == 'I' || p[0] == 'i') &&
            (p[1] == 'N' || p[1] == 'n') &&
            (p[2] == 'F' || p[2] == 'f')) {
            values[count++] = INF;
            p += 3;
        } else if (*p == '-' || isdigit((unsigned char)*p)) {
            values[count++] = atoi(p);
            if (*p == '-') p++;
            while (isdigit((unsigned char)*p)) p++;
        } else {
            // Skip unknown characters
            p++;
        }
    }
    return count;
}

// Parse a 7x7 matrix (for stack, etc.)
static void parse_7x7_matrix(FILE* fp, int matrix[NBPAIRS+1][NBPAIRS+1]) {
    char line[1024];
    int values[8];
    for (int i = 1; i <= NBPAIRS; i++) {
        if (!read_line(fp, line, sizeof(line))) break;
        int n = parse_int_line(line, values, 7);
        for (int j = 0; j < n && j < NBPAIRS; j++) {
            matrix[i][j+1] = values[j];
        }
    }
}

// Parse a 7x5x5 mismatch matrix
static void parse_mismatch_matrix(FILE* fp, int matrix[NBPAIRS+1][5][5]) {
    char line[1024];
    int values[8];
    // 7 base pair types * 5 rows each = 35 lines
    for (int bp = 1; bp <= NBPAIRS; bp++) {
        for (int i = 0; i < 5; i++) {
            if (!read_line(fp, line, sizeof(line))) break;
            int n = parse_int_line(line, values, 5);
            for (int j = 0; j < n && j < 5; j++) {
                matrix[bp][i][j] = values[j];
            }
        }
    }
}

// Parse dangle arrays (7x5)
static void parse_dangle(FILE* fp, int dangle[NBPAIRS+1][5]) {
    char line[1024];
    int values[8];
    for (int bp = 1; bp <= NBPAIRS; bp++) {
        if (!read_line(fp, line, sizeof(line))) break;
        int n = parse_int_line(line, values, 5);
        for (int j = 0; j < n && j < 5; j++) {
            dangle[bp][j] = values[j];
        }
    }
}

// Parse int11 (7x7x5x5)
static void parse_int11(FILE* fp, int int11[NBPAIRS+1][NBPAIRS+1][5][5]) {
    char line[1024];
    int values[8];
    for (int i = 1; i <= NBPAIRS; i++) {
        for (int j = 1; j <= NBPAIRS; j++) {
            for (int k = 0; k < 5; k++) {
                if (!read_line(fp, line, sizeof(line))) return;
                int n = parse_int_line(line, values, 5);
                for (int l = 0; l < n && l < 5; l++) {
                    int11[i][j][k][l] = values[l];
                }
            }
        }
    }
}

// Parse int21 (7x7x5x5x5)
static void parse_int21(FILE* fp, int int21[NBPAIRS+1][NBPAIRS+1][5][5][5]) {
    char line[1024];
    int values[8];
    for (int i = 1; i <= NBPAIRS; i++) {
        for (int j = 1; j <= NBPAIRS; j++) {
            for (int k = 0; k < 5; k++) {
                for (int l = 0; l < 5; l++) {
                    if (!read_line(fp, line, sizeof(line))) return;
                    int n = parse_int_line(line, values, 5);
                    for (int m = 0; m < n && m < 5; m++) {
                        int21[i][j][k][l][m] = values[m];
                    }
                }
            }
        }
    }
}

// Parse int22 (7x7x5x5x5x5) - but stored differently in file
static void parse_int22(FILE* fp, int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5]) {
    char line[1024];
    int values[8];
    // int22 is stored as blocks for each (i,j) pair
    for (int i = 1; i <= NBPAIRS; i++) {
        for (int j = 1; j <= NBPAIRS; j++) {
            for (int k = 0; k < 5; k++) {
                for (int l = 0; l < 5; l++) {
                    for (int m = 0; m < 5; m++) {
                        if (!read_line(fp, line, sizeof(line))) return;
                        int n = parse_int_line(line, values, 5);
                        for (int o = 0; o < n && o < 5; o++) {
                            int22[i][j][k][l][m][o] = values[o];
                        }
                    }
                }
            }
        }
    }
}

// Parse hairpin/bulge/interior arrays (size 31)
static void parse_loop_array(FILE* fp, int* arr, int size) {
    char line[1024];
    int values[32];
    int total = 0;
    while (total < size) {
        if (!read_line(fp, line, sizeof(line))) break;
        int n = parse_int_line(line, values, size - total);
        for (int i = 0; i < n && total < size; i++) {
            arr[total++] = values[i];
        }
    }
}

// Parse ML_params
static void parse_ml_params(FILE* fp) {
    char line[1024];
    int values[16];

    // First line: ML unpaired, ML closing, ML intern
    if (read_line(fp, line, sizeof(line))) {
        int n = parse_int_line(line, values, 8);
        if (n >= 1) raw_params.MLbase = values[0];
        if (n >= 2) raw_params.MLclosing = values[1];
        if (n >= 3) {
            for (int i = 0; i <= NBPAIRS; i++) {
                raw_params.MLintern[i] = values[2];
            }
        }
    }

    // Second line: enthalpies
    if (read_line(fp, line, sizeof(line))) {
        int n = parse_int_line(line, values, 8);
        if (n >= 1) raw_params.MLbase_enthalpies = values[0];
        if (n >= 2) raw_params.MLclosing_enthalpies = values[1];
        if (n >= 3) {
            for (int i = 0; i <= NBPAIRS; i++) {
                raw_params.MLintern_enthalpies[i] = values[2];
            }
        }
    }
}

// Parse NINIO
static void parse_ninio(FILE* fp) {
    char line[1024];
    int values[8];

    if (read_line(fp, line, sizeof(line))) {
        int n = parse_int_line(line, values, 5);
        for (int i = 0; i < n && i < 5; i++) {
            raw_params.ninio[i] = values[i];
        }
    }
    if (read_line(fp, line, sizeof(line))) {
        int n = parse_int_line(line, values, 5);
        for (int i = 0; i < n && i < 5; i++) {
            raw_params.ninio_enthalpies[i] = values[i];
        }
    }
}

// Parse Misc
static void parse_misc(FILE* fp) {
    char line[1024];
    int values[8];

    // First line: DuplexInit, TerminalAU, lxc
    if (read_line(fp, line, sizeof(line))) {
        int n = parse_int_line(line, values, 3);
        if (n >= 1) raw_params.DuplexInit = values[0];
        if (n >= 2) raw_params.TerminalAU = values[1];
        // lxc is a float, need special handling
        char* p = line;
        int count = 0;
        while (*p && count < 3) {
            while (*p && (isspace((unsigned char)*p) || *p == ',')) p++;
            if (*p == '\0' || *p == '/' || *p == '#') break;
            if (*p == '-' || *p == '.' || isdigit((unsigned char)*p)) {
                if (count == 2) {
                    raw_params.lxc = atof(p);
                }
                count++;
                while (*p && (*p == '-' || *p == '.' || isdigit((unsigned char)*p))) p++;
            } else {
                p++;
            }
        }
    }

    // Second line: enthalpies
    if (read_line(fp, line, sizeof(line))) {
        int n = parse_int_line(line, values, 2);
        if (n >= 1) raw_params.DuplexInit_enthalpies = values[0];
        if (n >= 2) raw_params.TerminalAU_enthalpies = values[1];
    }
}

// Parse special loops (Triloops, Tetraloops, Hexaloops)
static void parse_special_loops(FILE* fp, char* loops, int* energies, int* enthalpies,
                                 int loop_size, int max_loops, int* count) {
    char line[1024];
    *count = 0;

    while (read_line(fp, line, sizeof(line))) {
        char* p = skip_whitespace(line);
        if (*p == '#') break; // Next section

        // Format: SEQUENCE ENERGY [ENTHALPY]
        char seq[16];
        int energy = 0, enthalpy = 0;
        int n = sscanf(line, "%15s %d %d", seq, &energy, &enthalpy);
        if (n >= 2 && *count < max_loops) {
            strncpy(loops + (*count) * (loop_size + 1), seq, loop_size + 1);
            energies[*count] = energy;
            if (n >= 3) enthalpies[*count] = enthalpy;
            (*count)++;
        }
    }
}

void read_parameter_file(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open parameter file: %s\n", filename);
        exit(1);
    }

    // Initialize to zero
    memset(&raw_params, 0, sizeof(raw_params));
    raw_params.lxc = 107.856;  // Default value

    char line[1024];
    char section[256] = "";

    while (fgets(line, sizeof(line), fp)) {
        char* p = skip_whitespace(line);

        // Check for section header
        if (*p == '#' && p[1] != '#') {
            sscanf(p + 1, "%255s", section);

            if (strcmp(section, "stack") == 0) {
                parse_7x7_matrix(fp, raw_params.stack);
            } else if (strcmp(section, "stack_enthalpies") == 0) {
                parse_7x7_matrix(fp, raw_params.stack_enthalpies);
            } else if (strcmp(section, "mismatch_hairpin") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatchH);
            } else if (strcmp(section, "mismatch_hairpin_enthalpies") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatchH_enthalpies);
            } else if (strcmp(section, "mismatch_interior") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatchI);
            } else if (strcmp(section, "mismatch_interior_enthalpies") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatchI_enthalpies);
            } else if (strcmp(section, "mismatch_interior_1n") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatch1nI);
            } else if (strcmp(section, "mismatch_interior_1n_enthalpies") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatch1nI_enthalpies);
            } else if (strcmp(section, "mismatch_interior_23") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatch23I);
            } else if (strcmp(section, "mismatch_interior_23_enthalpies") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatch23I_enthalpies);
            } else if (strcmp(section, "mismatch_multi") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatchM);
            } else if (strcmp(section, "mismatch_multi_enthalpies") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatchM_enthalpies);
            } else if (strcmp(section, "mismatch_exterior") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatchExt);
            } else if (strcmp(section, "mismatch_exterior_enthalpies") == 0) {
                parse_mismatch_matrix(fp, raw_params.mismatchExt_enthalpies);
            } else if (strcmp(section, "dangle5") == 0) {
                parse_dangle(fp, raw_params.dangle5);
            } else if (strcmp(section, "dangle5_enthalpies") == 0) {
                parse_dangle(fp, raw_params.dangle5_enthalpies);
            } else if (strcmp(section, "dangle3") == 0) {
                parse_dangle(fp, raw_params.dangle3);
            } else if (strcmp(section, "dangle3_enthalpies") == 0) {
                parse_dangle(fp, raw_params.dangle3_enthalpies);
            } else if (strcmp(section, "int11") == 0) {
                parse_int11(fp, raw_params.int11);
            } else if (strcmp(section, "int11_enthalpies") == 0) {
                parse_int11(fp, raw_params.int11_enthalpies);
            } else if (strcmp(section, "int21") == 0) {
                parse_int21(fp, raw_params.int21);
            } else if (strcmp(section, "int21_enthalpies") == 0) {
                parse_int21(fp, raw_params.int21_enthalpies);
            } else if (strcmp(section, "int22") == 0) {
                parse_int22(fp, raw_params.int22);
            } else if (strcmp(section, "int22_enthalpies") == 0) {
                parse_int22(fp, raw_params.int22_enthalpies);
            } else if (strcmp(section, "hairpin") == 0) {
                parse_loop_array(fp, raw_params.hairpin, 31);
            } else if (strcmp(section, "hairpin_enthalpies") == 0) {
                parse_loop_array(fp, raw_params.hairpin_enthalpies, 31);
            } else if (strcmp(section, "bulge") == 0) {
                parse_loop_array(fp, raw_params.bulge, MAXLOOP + 1);
            } else if (strcmp(section, "bulge_enthalpies") == 0) {
                parse_loop_array(fp, raw_params.bulge_enthalpies, MAXLOOP + 1);
            } else if (strcmp(section, "interior") == 0) {
                parse_loop_array(fp, raw_params.internal_loop, MAXLOOP + 1);
            } else if (strcmp(section, "interior_enthalpies") == 0) {
                parse_loop_array(fp, raw_params.internal_loop_enthalpies, MAXLOOP + 1);
            } else if (strcmp(section, "ML_params") == 0) {
                parse_ml_params(fp);
            } else if (strcmp(section, "NINIO") == 0) {
                parse_ninio(fp);
            } else if (strcmp(section, "Misc") == 0) {
                parse_misc(fp);
            } else if (strcmp(section, "Triloops") == 0) {
                // Special loops need different parsing - they have sequence + energy format
                raw_params.num_triloops = 0;
                while (fgets(line, sizeof(line), fp)) {
                    p = skip_whitespace(line);
                    if (*p == '#') break;
                    char seq[16];
                    int energy = 0;
                    if (sscanf(line, "%15s %d", seq, &energy) >= 2) {
                        if (raw_params.num_triloops < 40) {
                            strncpy(raw_params.Triloops + raw_params.num_triloops * 6, seq, 5);
                            raw_params.Triloops[raw_params.num_triloops * 6 + 5] = '\0';
                            raw_params.Triloop_E[raw_params.num_triloops] = energy;
                            raw_params.num_triloops++;
                        }
                    }
                }
            } else if (strcmp(section, "Tetraloops") == 0) {
                raw_params.num_tetraloops = 0;
                while (fgets(line, sizeof(line), fp)) {
                    p = skip_whitespace(line);
                    if (*p == '#') break;
                    char seq[16];
                    int energy = 0;
                    if (sscanf(line, "%15s %d", seq, &energy) >= 2) {
                        if (raw_params.num_tetraloops < 200) {
                            strncpy(raw_params.Tetraloops + raw_params.num_tetraloops * 7, seq, 6);
                            raw_params.Tetraloops[raw_params.num_tetraloops * 7 + 6] = '\0';
                            raw_params.Tetraloop_E[raw_params.num_tetraloops] = energy;
                            raw_params.num_tetraloops++;
                        }
                    }
                }
            } else if (strcmp(section, "Hexaloops") == 0) {
                raw_params.num_hexaloops = 0;
                while (fgets(line, sizeof(line), fp)) {
                    p = skip_whitespace(line);
                    if (*p == '#') break;
                    char seq[16];
                    int energy = 0;
                    if (sscanf(line, "%15s %d", seq, &energy) >= 2) {
                        if (raw_params.num_hexaloops < 40) {
                            strncpy(raw_params.Hexaloops + raw_params.num_hexaloops * 9, seq, 8);
                            raw_params.Hexaloops[raw_params.num_hexaloops * 9 + 8] = '\0';
                            raw_params.Hexaloop_E[raw_params.num_hexaloops] = energy;
                            raw_params.num_hexaloops++;
                        }
                    }
                }
            } else if (strcmp(section, "END") == 0) {
                break;
            }
        }
    }

    fclose(fp);
    raw_params.loaded = 1;
}

// Temperature scaling: E(T) = H - T*S, where S = (H - E(37)) / 310.15
static int scale_energy(int energy37, int enthalpy, double temp_k) {
    if (energy37 == INF || enthalpy == INF) return INF;
    if (energy37 == FORBIDDEN || enthalpy == FORBIDDEN) return FORBIDDEN;

    double TT = (temp_k + K0) / (37.0 + K0);
    return (int)(enthalpy - (enthalpy - energy37) * TT + 0.5);
}

paramT* scale_parameters(void) {
    if (!raw_params.loaded) {
        fprintf(stderr, "Error: Parameters not loaded. Call read_parameter_file first.\n");
        exit(1);
    }

    paramT* P = (paramT*)calloc(1, sizeof(paramT));
    if (!P) {
        fprintf(stderr, "Error: Cannot allocate memory for parameters.\n");
        exit(1);
    }

    double temp_k = temperature; // Temperature in Celsius
    double TT = (temp_k + K0) / (37.0 + K0);

    P->temperature = temperature;
    P->lxc = raw_params.lxc * TT;

    // Scale stack energies
    for (int i = 0; i <= NBPAIRS; i++) {
        for (int j = 0; j <= NBPAIRS; j++) {
            P->stack[i][j] = scale_energy(raw_params.stack[i][j],
                                          raw_params.stack_enthalpies[i][j], temp_k);
        }
    }

    // Scale hairpin
    for (int i = 0; i < 31; i++) {
        P->hairpin[i] = scale_energy(raw_params.hairpin[i],
                                     raw_params.hairpin_enthalpies[i], temp_k);
    }

    // Scale bulge
    for (int i = 0; i <= MAXLOOP; i++) {
        P->bulge[i] = scale_energy(raw_params.bulge[i],
                                   raw_params.bulge_enthalpies[i], temp_k);
    }

    // Scale internal loop
    for (int i = 0; i <= MAXLOOP; i++) {
        P->internal_loop[i] = scale_energy(raw_params.internal_loop[i],
                                           raw_params.internal_loop_enthalpies[i], temp_k);
    }

    // Scale mismatch matrices
    for (int i = 0; i <= NBPAIRS; i++) {
        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 5; k++) {
                P->mismatchH[i][j][k] = scale_energy(raw_params.mismatchH[i][j][k],
                                                     raw_params.mismatchH_enthalpies[i][j][k], temp_k);
                P->mismatchI[i][j][k] = scale_energy(raw_params.mismatchI[i][j][k],
                                                     raw_params.mismatchI_enthalpies[i][j][k], temp_k);
                P->mismatch1nI[i][j][k] = scale_energy(raw_params.mismatch1nI[i][j][k],
                                                       raw_params.mismatch1nI_enthalpies[i][j][k], temp_k);
                P->mismatch23I[i][j][k] = scale_energy(raw_params.mismatch23I[i][j][k],
                                                       raw_params.mismatch23I_enthalpies[i][j][k], temp_k);
                P->mismatchM[i][j][k] = scale_energy(raw_params.mismatchM[i][j][k],
                                                     raw_params.mismatchM_enthalpies[i][j][k], temp_k);
                P->mismatchExt[i][j][k] = scale_energy(raw_params.mismatchExt[i][j][k],
                                                       raw_params.mismatchExt_enthalpies[i][j][k], temp_k);
            }
        }
    }

    // Scale dangles
    for (int i = 0; i <= NBPAIRS; i++) {
        for (int j = 0; j < 5; j++) {
            P->dangle5[i][j] = scale_energy(raw_params.dangle5[i][j],
                                            raw_params.dangle5_enthalpies[i][j], temp_k);
            P->dangle3[i][j] = scale_energy(raw_params.dangle3[i][j],
                                            raw_params.dangle3_enthalpies[i][j], temp_k);
        }
    }

    // Scale int11
    for (int i = 0; i <= NBPAIRS; i++) {
        for (int j = 0; j <= NBPAIRS; j++) {
            for (int k = 0; k < 5; k++) {
                for (int l = 0; l < 5; l++) {
                    P->int11[i][j][k][l] = scale_energy(raw_params.int11[i][j][k][l],
                                                        raw_params.int11_enthalpies[i][j][k][l], temp_k);
                }
            }
        }
    }

    // Scale int21
    for (int i = 0; i <= NBPAIRS; i++) {
        for (int j = 0; j <= NBPAIRS; j++) {
            for (int k = 0; k < 5; k++) {
                for (int l = 0; l < 5; l++) {
                    for (int m = 0; m < 5; m++) {
                        P->int21[i][j][k][l][m] = scale_energy(raw_params.int21[i][j][k][l][m],
                                                               raw_params.int21_enthalpies[i][j][k][l][m], temp_k);
                    }
                }
            }
        }
    }

    // Scale int22
    for (int i = 0; i <= NBPAIRS; i++) {
        for (int j = 0; j <= NBPAIRS; j++) {
            for (int k = 0; k < 5; k++) {
                for (int l = 0; l < 5; l++) {
                    for (int m = 0; m < 5; m++) {
                        for (int n = 0; n < 5; n++) {
                            P->int22[i][j][k][l][m][n] = scale_energy(raw_params.int22[i][j][k][l][m][n],
                                                                      raw_params.int22_enthalpies[i][j][k][l][m][n], temp_k);
                        }
                    }
                }
            }
        }
    }

    // Scale ninio
    for (int i = 0; i < 5; i++) {
        P->ninio[i] = scale_energy(raw_params.ninio[i], raw_params.ninio_enthalpies[i], temp_k);
    }

    // Scale ML parameters
    P->MLbase = scale_energy(raw_params.MLbase, raw_params.MLbase_enthalpies, temp_k);
    P->MLclosing = scale_energy(raw_params.MLclosing, raw_params.MLclosing_enthalpies, temp_k);
    for (int i = 0; i <= NBPAIRS; i++) {
        P->MLintern[i] = scale_energy(raw_params.MLintern[i], raw_params.MLintern_enthalpies[i], temp_k);
    }

    // Scale misc
    P->TerminalAU = scale_energy(raw_params.TerminalAU, raw_params.TerminalAU_enthalpies, temp_k);
    P->DuplexInit = scale_energy(raw_params.DuplexInit, raw_params.DuplexInit_enthalpies, temp_k);

    // Copy special loops (these are not temperature-scaled in the same way)
    memcpy(P->Tetraloops, raw_params.Tetraloops, sizeof(P->Tetraloops));
    memcpy(P->Tetraloop_E, raw_params.Tetraloop_E, sizeof(P->Tetraloop_E));
    memcpy(P->Triloops, raw_params.Triloops, sizeof(P->Triloops));
    memcpy(P->Triloop_E, raw_params.Triloop_E, sizeof(P->Triloop_E));
    memcpy(P->Hexaloops, raw_params.Hexaloops, sizeof(P->Hexaloops));
    memcpy(P->Hexaloop_E, raw_params.Hexaloop_E, sizeof(P->Hexaloop_E));

    // Initialize model details
    set_model_details(&P->model_details);

    return P;
}

void set_model_details(model_detailsT* md) {
    md->dangles = 2;
    md->special_hp = 1;
    md->noLP = 0;
    md->noGU = 0;
    md->noGUclosure = 0;
    md->logML = 0;
    md->circ = 0;
    md->gquad = 0;
}
