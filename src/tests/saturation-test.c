/* libh2o -- steam & water properties
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
	CHECK(steam_saturation_p_T(300), 0.353658941E-2, 1E-11);
	CHECK(steam_saturation_p_T(500), 0.263889776E+1, 1E-8);
	CHECK(steam_saturation_p_T(600), 0.123443146E+2, 1E-7);

	return 0;
}
