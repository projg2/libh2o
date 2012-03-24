/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region1.h"
#include "region2.h"
#include "saturation.h"

/* Create values by interpolating R1 & R2
 * for p <= 10 MPa, use R2-meta instead */

double h2o_region4_v_Tx(double T, double x) /* [K, 0..1] -> [m³/kg] */
{
	double p = h2o_saturation_p_T(T);

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
	double p = h2o_saturation_p_T(T);

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
	double p = h2o_saturation_p_T(T);

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
	double p = h2o_saturation_p_T(T);

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
	double p = h2o_saturation_p_T(T);

	double s1 = h2o_region1_s_pT(p, T);
	double s2 = p <= 10 ? h2o_region2_meta_s_pT(p, T)
			: h2o_region2_s_pT(p, T);

	return (s - s1) / (s2 - s1);
}

double h2o_region4_x_Th(double T, double h) /* [K, kJ/kg] -> [0..1] */
{
	double p = h2o_saturation_p_T(T);

	double h1 = h2o_region1_h_pT(p, T);
	double h2 = p <= 10 ? h2o_region2_meta_h_pT(p, T)
			: h2o_region2_h_pT(p, T);

	return (h - h1) / (h2 - h1);
}
