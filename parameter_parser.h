#ifndef PARAMETER_PARSER_H
#define PARAMETER_PARSER_H

#include "params.h"
#include "memory_types.h"
#include <stdexcept>

/**
 * Read energy parameters from a ViennaRNA-format parameter file.
 * This replaces the ViennaRNA library dependency.
 *
 * @param filename Path to the parameter file (e.g., rna_turner2004.par)
 * @throws std::runtime_error if file cannot be opened or parsed
 */
void read_parameter_file(const char* filename);

/**
 * Scale the energy parameters for the specified temperature.
 *
 * @param temperature Temperature in Celsius
 * @return Smart pointer to newly allocated paramT structure with scaled parameters
 * @throws std::runtime_error if parameters not loaded or allocation fails
 */
fftbor::ParamPtr scale_parameters(double temperature);

/**
 * Initialize model details with default values.
 *
 * @param md Pointer to model_detailsT structure to initialize
 */
void set_model_details(model_detailsT* md);

/**
 * Load parameters and scale them for the given context.
 * Convenience function that calls read_parameter_file and scale_parameters.
 *
 * @param ctx Context with energy_file and temperature set
 * @throws std::runtime_error on any error
 */
void load_parameters(fftbor::Context& ctx);

#endif /* PARAMETER_PARSER_H */
