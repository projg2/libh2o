/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <stdlib.h>

#include "region2.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 6.3: Backward Equations */

static const double n[] = {
	+0.00000000000000E+0,

	+0.90584278514723E+3,
	-0.67955786399241E+0,
	+0.12809002730136E-3,
	+0.26526571908428E+4,
	+0.45257578905948E+1
};

double h2o_region2_b2bc_p_h(double h) /* [kJ/kg] -> [MPa] */
{
	return quadr_value(n[3], n[2], n[1], h);
}

double h2o_region2_b2bc_h_p(double p) /* [MPa] -> [kJ/kg] */
{
	return sqrt((p - n[5]) / n[3]) + n[4];
}

enum h2o_region2_subregion h2o_region2_subregion_ph(double p, double h)
{
	if (p < 4)
		return H2O_REGION2A;
	if (p < h2o_region2_b2bc_p_h(h))
		return H2O_REGION2B;
	else
		return H2O_REGION2C;
}

enum h2o_region2_subregion h2o_region2_subregion_ps(double p, double s)
{
	if (p < 4)
		return H2O_REGION2A;
	if (s >= 5.85)
		return H2O_REGION2B;
	else
		return H2O_REGION2C;
}

typedef double (*twoarg_func_t)(double, double);

/* this function intends to make compiler happy. */
static double impossible_happened(double a, double b)
{
	abort();

	if (a || b)
	{
	}
}

double h2o_region2_T_ph(double p, double h)
{
	twoarg_func_t T_getter;

	switch (h2o_region2_subregion_ph(p, h))
	{
		case H2O_REGION2A:
			T_getter = &h2o_region2a_T_ph;
			break;
		case H2O_REGION2B:
			T_getter = &h2o_region2b_T_ph;
			break;
		case H2O_REGION2C:
			T_getter = &h2o_region2c_T_ph;
			break;
		default:
			T_getter = &impossible_happened;
	}

	return T_getter(p, h);
}

double h2o_region2_T_ps(double p, double s)
{
	twoarg_func_t T_getter;

	switch (h2o_region2_subregion_ps(p, s))
	{
		case H2O_REGION2A:
			T_getter = &h2o_region2a_T_ps;
			break;
		case H2O_REGION2B:
			T_getter = &h2o_region2b_T_ps;
			break;
		case H2O_REGION2C:
			T_getter = &h2o_region2c_T_ps;
			break;
		default:
			T_getter = &impossible_happened;
	}

	return T_getter(p, s);
}