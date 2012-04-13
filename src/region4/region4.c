/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "consts.h"
#include "region1.h"
#include "region2.h"
#include "region3.h"
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

typedef double (*twoarg_func_t)(double, double);

/* Create values by interpolating R1 & R2
 * for p <= 10 MPa, use R2-meta instead */

static inline double region4_interp(
		twoarg_func_t sat_water_func,
		twoarg_func_t sat_steam_func,
		double T, double x)
{
	double p = h2o_region4_p_T(T);

	if (x == 0)
		return sat_water_func(p, T);
	else
	{
		double v2 = sat_steam_func(p, T);

		if (x == 1)
			return v2;
		else
		{
			double v1 = sat_water_func(p, T);

			return v1 + (v2 - v1) * x;
		}
	}
}

enum h2o_region4_subregion {
	H2O_REGION4_12META,
	H2O_REGION4_12,
	H2O_REGION4_3CT,
	H2O_REGION4_3ST,
	H2O_REGION4_3RS,
	H2O_REGION4_3UX,
	H2O_REGION4_3UZ,
	H2O_REGION4_3YZ
};

static enum h2o_region4_subregion h2o_region4_subregion_T(double T)
{
	if (T <= Tsat2metamax)
		return H2O_REGION4_12META;
	else if (T <= Tsat12max)
		return H2O_REGION4_12;
	else if (T <= Tsat3cmax)
		return H2O_REGION4_3CT;
	else if (T <= Tsat3tmax)
		return H2O_REGION4_3ST;
	else if (T <= Tsat3rsmax)
		return H2O_REGION4_3RS;
	else if (T <= Tsat3xmax)
		return H2O_REGION4_3UX;
	else if (T <= Tsat3umax)
		return H2O_REGION4_3UZ;
	else
		return H2O_REGION4_3YZ;
}

double h2o_region4_v_Tx(double T, double x) /* [K, 0..1] -> [m³/kg] */
{
	twoarg_func_t water_func, steam_func;

	switch (h2o_region4_subregion_T(T))
	{
		case H2O_REGION4_12META:
			water_func = h2o_region1_v_pT;
			steam_func = h2o_region2_meta_v_pT;
			break;
		case H2O_REGION4_12:
			water_func = h2o_region1_v_pT;
			steam_func = h2o_region2_v_pT;
			break;
		case H2O_REGION4_3CT:
			water_func = h2o_region3c_v_pT;
			steam_func = h2o_region3t_v_pT;
			break;
		case H2O_REGION4_3ST:
			water_func = h2o_region3s_v_pT;
			steam_func = h2o_region3t_v_pT;
			break;
		case H2O_REGION4_3RS:
			water_func = h2o_region3s_v_pT;
			steam_func = h2o_region3r_v_pT;
			break;
		case H2O_REGION4_3UX:
			water_func = h2o_region3u_v_pT;
			steam_func = h2o_region3x_v_pT;
			break;
		case H2O_REGION4_3UZ:
			water_func = h2o_region3u_v_pT;
			steam_func = h2o_region3z_v_pT;
			break;
		case H2O_REGION4_3YZ:
			water_func = h2o_region3y_v_pT;
			steam_func = h2o_region3z_v_pT;
			break;
	}

	return region4_interp(water_func, steam_func, T, x);
}

static double sat_region3_u1_pT(double p, double T)
{
	double v = h2o_region4_v_Tx(T, 0);

	return h2o_region3_u_rhoT(1/v, T);
}

static double sat_region3_u2_pT(double p, double T)
{
	double v = h2o_region4_v_Tx(T, 1);

	return h2o_region3_u_rhoT(1/v, T);
}

double h2o_region4_u_Tx(double T, double x) /* [K, 0..1] -> [kJ/kg] */
{
	twoarg_func_t water_func, steam_func;

	switch (h2o_region4_subregion_T(T))
	{
		case H2O_REGION4_12META:
			water_func = h2o_region1_u_pT;
			steam_func = h2o_region2_meta_u_pT;
			break;
		case H2O_REGION4_12:
			water_func = h2o_region1_u_pT;
			steam_func = h2o_region2_u_pT;
			break;
		default:
			water_func = sat_region3_u1_pT;
			steam_func = sat_region3_u2_pT;
	}

	return region4_interp(water_func, steam_func, T, x);
}

static double sat_region3_s1_pT(double p, double T)
{
	double v = h2o_region4_v_Tx(T, 0);

	return h2o_region3_s_rhoT(1/v, T);
}

static double sat_region3_s2_pT(double p, double T)
{
	double v = h2o_region4_v_Tx(T, 1);

	return h2o_region3_s_rhoT(1/v, T);
}

double h2o_region4_s_Tx(double T, double x) /* [K, 0..1] -> [kJ/kgK] */
{
	twoarg_func_t water_func, steam_func;

	switch (h2o_region4_subregion_T(T))
	{
		case H2O_REGION4_12META:
			water_func = h2o_region1_s_pT;
			steam_func = h2o_region2_meta_s_pT;
			break;
		case H2O_REGION4_12:
			water_func = h2o_region1_s_pT;
			steam_func = h2o_region2_s_pT;
			break;
		default:
			water_func = sat_region3_s1_pT;
			steam_func = sat_region3_s2_pT;
	}

	return region4_interp(water_func, steam_func, T, x);
}

static double sat_region3_h1_pT(double p, double T)
{
	double v = h2o_region4_v_Tx(T, 0);

	return h2o_region3_h_rhoT(1/v, T);
}

static double sat_region3_h2_pT(double p, double T)
{
	double v = h2o_region4_v_Tx(T, 1);

	return h2o_region3_h_rhoT(1/v, T);
}

double h2o_region4_h_Tx(double T, double x) /* [K, 0..1] -> [kJ/kg] */
{
	twoarg_func_t water_func, steam_func;

	switch (h2o_region4_subregion_T(T))
	{
		case H2O_REGION4_12META:
			water_func = h2o_region1_h_pT;
			steam_func = h2o_region2_meta_h_pT;
			break;
		case H2O_REGION4_12:
			water_func = h2o_region1_h_pT;
			steam_func = h2o_region2_h_pT;
			break;
		default:
			water_func = sat_region3_h1_pT;
			steam_func = sat_region3_h2_pT;
	}

	return region4_interp(water_func, steam_func, T, x);
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
