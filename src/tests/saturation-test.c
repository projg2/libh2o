/* libh2o -- h2o & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "tests.h"
#include "saturation.h"

int main(void)
{
	CHECK(h2o_saturation_p_T(300), 0.353658941E-2, 1E-11);
	CHECK(h2o_saturation_p_T(500), 0.263889776E+1, 1E-8);
	CHECK(h2o_saturation_p_T(600), 0.123443146E+2, 1E-7);

	CHECK(h2o_saturation_T_p(0.1), 0.372755919E3, 1E-6);
	CHECK(h2o_saturation_T_p(1.0), 0.453035632E3, 1E-6);
	CHECK(h2o_saturation_T_p(10.), 0.584149488E3, 1E-6);

	return 0;
}
