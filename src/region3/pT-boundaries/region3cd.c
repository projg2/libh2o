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
	+0.000000000000000E+0,

	+0.585276966696349E+3, +0.278233532206915E+1,
	-0.127283549295878E-1, +0.159090746562729E-3
};

double h2o_region3cd_T_p(double p)
{
	return poly_value(p, 0, 3, 0, n);
}
