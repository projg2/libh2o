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

	+0.918419702359447E3, -0.191887498864292E4,
	+0.154793642129415E4, -0.187661219490113E3,
	+0.213144632222113E2
};

double h2o_region3ab_T_p(double p)
{
	return poly_value(log(p), -2, 3, 0, n);
}
