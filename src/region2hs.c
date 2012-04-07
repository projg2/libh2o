/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region2.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Pressure as a Function
 * of Enthalpy and Entropy p(h,s) to the IAPWS Industrial Formulation 1997
 * for the Thermodynamic Properties of Water and Steam;
 * s. 6: Backward Equation p(h,s) for Region 1 */

static const double n[] = {
	+0.000000000000000E0,

	-0.349898083432139E4,
	+0.257560716905876E4,
	-0.421073558227969E3,
	+0.276349063799944E2
};

double h2o_region2_b2ab_h_s(double s)
{
	return deg3_value(n[4], n[3], n[2], n[1], s);
}
