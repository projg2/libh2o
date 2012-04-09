/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"

enum h2o_region3_subregion h2o_region3_subregion_ph(double p, double h)
{
	if (h <= h2o_region3_b3ab_h_p(p))
		return H2O_REGION3A;
	else
		return H2O_REGION3B;
}

enum h2o_region3_subregion h2o_region3_subregion_ps(double p, double s)
{
	if (s <= 4.41202148223476) /* scrit */
		return H2O_REGION3A;
	else
		return H2O_REGION3B;
}
