/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <stdlib.h> /* abort() */

#include "region3.h"

enum h2o_region3_subregion h2o_region3_subregion_ph(double p, double h)
{
	if (h <= h2o_region3_b3ab_h_p(p))
		return H2O_REGION3A;
	else
		return H2O_REGION3B;
}

enum h2o_region3_subregion h2o_region3_subregion_ps(double p, double s)
{
	if (s <= 4.41202148223476) /* scrit */
		return H2O_REGION3A;
	else
		return H2O_REGION3B;
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

double h2o_region3_T_ph(double p, double h)
{
	twoarg_func_t T_getter;

	switch (h2o_region3_subregion_ph(p, h))
	{
		case H2O_REGION3A:
			T_getter = &h2o_region3a_T_ph;
			break;
		case H2O_REGION3B:
			T_getter = &h2o_region3b_T_ph;
			break;
		default:
			T_getter = &impossible_happened;
	}

	return T_getter(p, h);
}

double h2o_region3_v_ph(double p, double h)
{
	twoarg_func_t v_getter;

	switch (h2o_region3_subregion_ph(p, h))
	{
		case H2O_REGION3A:
			v_getter = &h2o_region3a_v_ph;
			break;
		case H2O_REGION3B:
			v_getter = &h2o_region3b_v_ph;
			break;
		default:
			v_getter = &impossible_happened;
	}

	return v_getter(p, h);
}

double h2o_region3_T_ps(double p, double s)
{
	twoarg_func_t T_getter;

	switch (h2o_region3_subregion_ps(p, s))
	{
		case H2O_REGION3A:
			T_getter = &h2o_region3a_T_ps;
			break;
		case H2O_REGION3B:
			T_getter = &h2o_region3b_T_ps;
			break;
		default:
			T_getter = &impossible_happened;
	}

	return T_getter(p, s);
}

double h2o_region3_v_ps(double p, double s)
{
	twoarg_func_t v_getter;

	switch (h2o_region3_subregion_ps(p, s))
	{
		case H2O_REGION3A:
			v_getter = &h2o_region3a_v_ps;
			break;
		case H2O_REGION3B:
			v_getter = &h2o_region3b_v_ps;
			break;
		default:
			v_getter = &impossible_happened;
	}

	return v_getter(p, s);
}
