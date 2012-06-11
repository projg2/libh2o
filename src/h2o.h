/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_H2O_H
#define _H2O_H2O_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

#include "boundaries.h"

/**
 * A struct describing IF97 state point.
 *
 * 1) Create a new one using h2o_new_*(),
 * 2) use h2o_is_valid() or manually ensure that
 * region != H2O_REGION_OUT_OF_RANGE,
 * 3) call any of h2o_get_*() on the structure.
 */
typedef struct
{
	enum h2o_region region;

/*private:*/
	double _arg1, _arg2;
} h2o_t;

/**
 * Property cheat sheet:
 *
 * p - pressure [MPa],
 * T - temperature [K],
 * u - specific internal energy [kJ/kg],
 * h - specific enthalpy [kJ/kg],
 * s - specific enthropy [kJ/kgK],
 * v - specific volume [m³/kg],
 * cp - specific isobaric heat capacity [kJ/kgK],
 * cv - specific isochoric heat capacity [kJ/kgK],
 * w - speed of sound [m/s],
 * rho - density [kg/m³],
 * x - dryness [0..1].
 */

/**
 * Initializers.
 *
 * Those functions create a new h2o_t (state point) using two given
 * parameters.
 *
 * If the parameters are out of range, the struct will have region set
 * to H2O_REGION_OUT_OF_RANGE. You can use h2o_is_valid() to easily
 * check for that. Such a struct must not be passed to h2o_get_*().
 */

h2o_t h2o_new_pT(double p, double T);
h2o_t h2o_new_Tx(double T, double x);
h2o_t h2o_new_px(double p, double x);
h2o_t h2o_new_ph(double p, double h);
h2o_t h2o_new_ps(double p, double s);
h2o_t h2o_new_hs(double h, double s);
h2o_t h2o_new_rhoT(double rho, double T);

/**
 * Check whether a particular state point is valid and in range.
 *
 * Returns a true (non-zero) value if it is, false (zero) otherwise.
 */

int h2o_is_valid(const h2o_t state);

/**
 * Get region for a particular state point.
 *
 * If the parameters were out of range, returns H2O_REGION_OUT_OF_RANGE.
 * Otherwise, returns one of the H2O_REGIONn.
 *
 * It is equivalent to reading the public state.region field.
 */

enum h2o_region h2o_get_region(const h2o_t state);

/**
 * Getters.
 *
 * Those functions get various properties for a given state point.
 * The state must be initialized (using h2o_new_*()) and valid (check
 * using h2o_is_valid()).
 *
 * Note: h2o_get_x() must not be used on Region 3 data. You should check
 * the region before using it.
 */

double h2o_get_p(const h2o_t state);
double h2o_get_T(const h2o_t state);
double h2o_get_x(const h2o_t state);
double h2o_get_rho(const h2o_t state);

double h2o_get_v(const h2o_t state);
double h2o_get_u(const h2o_t state);
double h2o_get_h(const h2o_t state);
double h2o_get_s(const h2o_t state);
double h2o_get_cp(const h2o_t state);
double h2o_get_cv(const h2o_t state);
double h2o_get_w(const h2o_t state);

/**
 * Perform an expansion calculation from the given state point.
 *
 * @pout: target pressure [MPa]
 * @eta: (optional) isenthropic efficiency [0..1]
 *
 * Returns a new state point. If the pre- or post-expansions parameters
 * were out of supported range, the struct will have region set
 * to H2O_REGION_OUT_OF_RNAGE.
 *
 * The variant without @eta assumes ideal expansion (@eta = 1).
 */

h2o_t h2o_expand(const h2o_t in_state, double pout);
h2o_t h2o_expand_real(const h2o_t in_state, double pout, double eta);

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_H2O_H*/
