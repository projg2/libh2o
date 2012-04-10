/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

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
