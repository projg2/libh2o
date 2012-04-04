/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region1.h"
#include "region2.h"
#include "region4.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 8: Equations for Region 4 */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
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
double h2o_region4_p_T(double T) /* p [MPa] = f(T [K]) */
{
	double theta = T + n[9] / (T - n[10]);

	double A = quadr_value(n[0], n[1], n[2], theta);
	double B = quadr_value(n[3], n[4], n[5], theta);
	double C = quadr_value(n[6], n[7], n[8], theta);

	/* gcc is not smart enough to notice 2*C being used twice */
	double twoC = 2 * C;

	double delta = pow2(B) - 2 * A * twoC;
	double ret = twoC / (-B + sqrt(delta));
	double retqu = pow4(ret);

	return retqu;
}

/* 8.2 The Saturation-Temperature Equation (Backward Equation) */
double h2o_region4_T_p(double p) /* T [K] = f(p [MPa]) */
{
	double beta = pow(p, 0.25);

	double E = quadr_value(n[0], n[3], n[6], beta);
	double F = quadr_value(n[1], n[4], n[7], beta);
	double G = quadr_value(n[2], n[5], n[8], beta);

	double delta = pow2(F) - 4 * E * G;
	double halfD = G / (-F - sqrt(delta));

	double subexpr = pow2(n[10] / 2) + pow2(halfD) - n[9] - n[10] * halfD;
	double ret = n[10] / 2 + halfD - sqrt(subexpr);

	return ret;
}

/* Create values by interpolating R1 & R2
 * for p <= 10 MPa, use R2-meta instead */

double h2o_region4_v_Tx(double T, double x) /* [K, 0..1] -> [m³/kg] */
{
	double p = h2o_region4_p_T(T);

	if (x == 0)
		return h2o_region1_v_pT(p, T);
	else
	{
		double v2 = p <= 10 ? h2o_region2_meta_v_pT(p, T)
				: h2o_region2_v_pT(p, T);

		if (x == 1)
			return v2;
		else
		{
			double v1 = h2o_region1_v_pT(p, T);

			return v1 + (v2 - v1) * x;
		}
	}
}

double h2o_region4_u_Tx(double T, double x) /* [K, 0..1] -> [kJ/kg] */
{
	double p = h2o_region4_p_T(T);

	if (x == 0)
		return h2o_region1_u_pT(p, T);
	else
	{
		double u2 = p <= 10 ? h2o_region2_meta_u_pT(p, T)
				: h2o_region2_u_pT(p, T);

		if (x == 1)
			return u2;
		else
		{
			double u1 = h2o_region1_u_pT(p, T);

			return u1 + (u2 - u1) * x;
		}
	}
}

double h2o_region4_s_Tx(double T, double x) /* [K, 0..1] -> [kJ/kgK] */
{
	double p = h2o_region4_p_T(T);

	if (x == 0)
		return h2o_region1_s_pT(p, T);
	else
	{
		double s2 = p <= 10 ? h2o_region2_meta_s_pT(p, T)
				: h2o_region2_s_pT(p, T);

		if (x == 1)
			return s2;
		else
		{
			double s1 = h2o_region1_s_pT(p, T);

			return s1 + (s2 - s1) * x;
		}
	}
}

double h2o_region4_h_Tx(double T, double x) /* [K, 0..1] -> [kJ/kg] */
{
	double p = h2o_region4_p_T(T);

	if (x == 0)
		return h2o_region1_h_pT(p, T);
	else
	{
		double h2 = p <= 10 ? h2o_region2_meta_h_pT(p, T)
				: h2o_region2_h_pT(p, T);

		if (x == 1)
			return h2;
		else
		{
			double h1 = h2o_region1_h_pT(p, T);

			return h1 + (h2 - h1) * x;
		}
	}
}

double h2o_region4_x_Ts(double T, double s) /* [K, kJ/kgK] -> [0..1] */
{
	double p = h2o_region4_p_T(T);

	double s1 = h2o_region1_s_pT(p, T);
	double s2 = p <= 10 ? h2o_region2_meta_s_pT(p, T)
			: h2o_region2_s_pT(p, T);

	return (s - s1) / (s2 - s1);
}

double h2o_region4_x_Th(double T, double h) /* [K, kJ/kg] -> [0..1] */
{
	double p = h2o_region4_p_T(T);

	double h1 = h2o_region1_h_pT(p, T);
	double h2 = p <= 10 ? h2o_region2_meta_h_pT(p, T)
			: h2o_region2_h_pT(p, T);

	return (h - h1) / (h2 - h1);
}
