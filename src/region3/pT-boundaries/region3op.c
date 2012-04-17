/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Specific Volume
 * as a Function of Pressure and Temperature v(p,T)
 * for Region 3 of the IAPWS Industrial Formulation 1997 for the
 * Thermodynamic Properties of Water and Steam */

/* reordered to sort powers incrementally */
static const double n[] = {
	+0.000000000000000E0,

	-0.152313732937084E4, +0.773845935768222E3,
	+0.969461372400213E3, -0.332500170441278E3,
	+0.642859598466067E2
};

double h2o_region3op_T_p(double p)
{
	return poly_value(log(p), -2, 2, 0, n);
}
