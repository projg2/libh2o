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
static const double n_store[] = {
	+0.918419702359447E3, -0.191887498864292E4,
	+0.154793642129415E4, -0.187661219490113E3,
	+0.213144632222113E2
};

static const double *n = &n_store[2];

double h2o_region3ab_T_p(double p)
{
	double pi = p;

	int i;
	double sum = 0;

	double pipowers_store[3 + 2];
		/* shift it for negative indices */
	double* pipowers = &pipowers_store[2];

	fill_powers_incr(pipowers, -2, 3, log(pi), 0);

	for (i = -2; i <= 2; ++i)
		sum += n[i] * pipowers[i];

	return sum;
}
