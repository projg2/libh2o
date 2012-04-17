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

	+0.584814781649163E+3, -0.616179320924617E+0,
	+0.260763050899562E+0, -0.587071076864459E-2,
	+0.515308185433082E-4
};

double h2o_region3ij_T_p(double p)
{
	return poly_value(p, 0, 4, 0, n);
}
