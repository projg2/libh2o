/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <math.h>

#include "saturation.h"

static double n[] = {
	+1.00000000000000E0,
	+0.11670521452767E4,
	-0.72421316703206E6,
	-0.17073846940092E2,
	+0.12020824702470E5,
	-0.32325550323333E7,
	+0.14915108613530E2,
	-0.48232657361591E4,
	+0.40511340542057E6,
	-0.23855557567849E0,
	+0.65017534844798E3
};

double h2o_saturation_p_T(double T) /* p [MPa] = f(T [K]) */
{
	double theta = T + n[9] / (T - n[10]);
	double thetasq = theta * theta;

	double A = n[0] * thetasq + n[1] * theta + n[2];
	double B = n[3] * thetasq + n[4] * theta + n[5];
	double C = n[6] * thetasq + n[7] * theta + n[8];

	double delta = B * B - 4 * A * C;
	double ret = 2 * C / (-B + sqrt(delta));
	double retqu = ret * ret * ret * ret;

	return retqu;
}
