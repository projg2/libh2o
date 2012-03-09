/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <math.h>

#include "boundaries.h"
#include "saturation.h"

/* Based on IF97-Rev, s. 2: Structure of the Formulation
 * and s. 4: Auxiliary Equation for the Boundary between Regions 2 and 3 */

static const double n[] = {
	+0.00000000000000E+0,

	+0.34805185628969E+3,
	-0.11671859879975E+1,
	+0.10192970039326E-2,
	+0.57254459862746E+3,
	+0.13918839778870E+2
};

double h2o_b23_p_T(double T) /* [K] -> [MPa] */
{
	return ((n[3] * T) + n[2]) * T + n[1];
}

double h2o_b23_T_p(double p) /* [MPa] -> [K] */
{
	return sqrt((p - n[5]) / n[3]) + n[4];
}

enum h2o_region h2o_region_pT(double p, double T) /* [MPa, K] */
{
	if (T < 273.15 || T > 2273.15 || p < 0 || p > 100)
		return H2O_REGION_OUT_OF_RANGE;

	else if (T < 623.15) /* 1 or 2 */
	{
		if (p >= h2o_saturation_p_T(T))
			return H2O_REGION1;
		else
			return H2O_REGION2;
	}

	else if (T < 1073.15) /* 3 or 2 */
	{
		if (p >= h2o_b23_p_T(T))
			return H2O_REGION3;
		else
			return H2O_REGION2;
	}

	else /* 5? */
	{
		if (p > 50)
			return H2O_REGION_OUT_OF_RANGE;
		else
			return H2O_REGION5;
	}
}
