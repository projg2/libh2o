/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "boundaries.h"
#include "xmath.h"

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
	return quadr_value(n[3], n[2], n[1], T);
}

double h2o_b23_T_p(double p) /* [MPa] -> [K] */
{
	return sqrt((p - n[5]) / n[3]) + n[4];
}
