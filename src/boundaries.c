/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "boundaries.h"
#include "saturation.h"

/* Based on IF97-Rev, s. 2: Structure of the Formulation */

enum h2o_region h2o_region_pT(double p, double T) /* [MPa, K] */
{
	if (T < 273.15 || T > 2273.15 || p < 0 || p > 100)
		return H2O_REGION_OUT_OF_RANGE;

	if (T < 623.15) /* 1 or 4 */
	{
		if (p >= h2o_saturation_p_T(T))
			return H2O_REGION1;
		else
			return H2O_REGION2;
	}

	if (T > 1073.15) /* maybe 5? */
	{
		if (p > 50)
			return H2O_REGION_OUT_OF_RANGE;
		else
			return H2O_REGION5;
	}

	/* XXX: 3/2 */

	return H2O_REGION_OUT_OF_RANGE;
}
