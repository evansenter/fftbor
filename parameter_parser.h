#ifndef PARAMETER_PARSER_H
#define PARAMETER_PARSER_H

#include "params.h"

/**
 * Read energy parameters from a ViennaRNA-format parameter file.
 * This replaces the ViennaRNA library dependency.
 *
 * @param filename Path to the parameter file (e.g., rna_turner2004.par)
 */
void read_parameter_file(const char* filename);

/**
 * Scale the energy parameters for the current temperature.
 * Temperature is expected to be set via the global 'temperature' variable.
 *
 * @return Pointer to newly allocated paramT structure with scaled parameters
 */
paramT* scale_parameters(void);

/**
 * Initialize model details with default values.
 *
 * @param md Pointer to model_detailsT structure to initialize
 */
void set_model_details(model_detailsT* md);

#endif /* PARAMETER_PARSER_H */
