/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Revised Supplementary Release on Backward Equations for the Functions
 * T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial
 * Formulation 1997 for the Thermodynamic Properties of Water and Steam
 * s. 3.2: Structure of the Equation Set */

static const double n[] = {
	+0.000000000000000E+0,

	+0.201464004206875E+4,
	+0.374696550136983E+1,
	-0.219921901054187E-1,
	+0.875131686009950E-4
};

double h2o_region3_b3ab_h_p(double p)
{
	return deg3_value(n[4], n[3], n[2], n[1], p);
}
