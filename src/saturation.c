/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <math.h>

#include "saturation.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 8: Equations for Region 4 */

/* coefficient table; n[0] added for convenience */
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

/* 8.1 The Saturation-Pressure Equation (Basic Equation) */
double h2o_saturation_p_T(double T) /* p [MPa] = f(T [K]) */
{
	double theta = T + n[9] / (T - n[10]);
	double thetasq = pow2(theta);

	double A = n[0] * thetasq + n[1] * theta + n[2];
	double B = n[3] * thetasq + n[4] * theta + n[5];
	double C = n[6] * thetasq + n[7] * theta + n[8];

	/* gcc is not smart enough to notice 2*C being used twice */
	double twoC = 2 * C;

	double delta = pow2(B) - 2 * A * twoC;
	double ret = twoC / (-B + sqrt(delta));
	double retqu = pow4(ret);

	return retqu;
}

/* 8.2 The Saturation-Temperature Equation (Backward Equation) */
double h2o_saturation_T_p(double p) /* T [K] = f(p [MPa]) */
{
	double betasq = sqrt(p);
	double beta = sqrt(betasq);

	double E = n[0] * betasq + n[3] * beta + n[6];
	double F = n[1] * betasq + n[4] * beta + n[7];
	double G = n[2] * betasq + n[5] * beta + n[8];

	/* gcc is not smart enough to notice 2*G being used twice */
	double twoG = 2 * G;

	double delta = pow2(F) - 2 * E * twoG;
	double D = twoG / (-F - sqrt(delta));

	double Dadj = n[10] + D;
	double subexpr = pow2(Dadj) - 4 * (n[9] + n[10] * D);

	double ret = (Dadj - sqrt(subexpr)) / 2;

	return ret;
}
