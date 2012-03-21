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

double h2o_region4_v_Tx(double T, double x) /* [K, 0..1] -> [m³/kg] */
{
	double p = h2o_saturation_p_T(T);

	if (x == 0)
		return h2o_region1_v_pT(p, T);
	else if (x == 1)
		return h2o_region2_v_pT(p, T);
	else
	{
		double v1 = h2o_region1_v_pT(p, T);
		double v2 = h2o_region2_v_pT(p, T);

		return v1 + (v2 - v1) * x;
	}
}

double h2o_region4_u_Tx(double T, double x) /* [K, 0..1] -> [kJ/kg] */
{
	double p = h2o_saturation_p_T(T);

	if (x == 0)
		return h2o_region1_u_pT(p, T);
	else if (x == 1)
		return h2o_region2_u_pT(p, T);
	else
	{
		double u1 = h2o_region1_u_pT(p, T);
		double u2 = h2o_region2_u_pT(p, T);

		return u1 + (u2 - u1) * x;
	}
}

double h2o_region4_s_Tx(double T, double x) /* [K, 0..1] -> [kJ/kgK] */
{
	double p = h2o_saturation_p_T(T);

	if (x == 0)
		return h2o_region1_s_pT(p, T);
	else if (x == 1)
		return h2o_region2_s_pT(p, T);
	else
	{
		double s1 = h2o_region1_s_pT(p, T);
		double s2 = h2o_region2_s_pT(p, T);

		return s1 + (s2 - s1) * x;
	}
}

double h2o_region4_h_Tx(double T, double x) /* [K, 0..1] -> [kJ/kg] */
{
	double p = h2o_saturation_p_T(T);

	if (x == 0)
		return h2o_region1_h_pT(p, T);
	else if (x == 1)
		return h2o_region2_h_pT(p, T);
	else
	{
		double h1 = h2o_region1_h_pT(p, T);
		double h2 = h2o_region2_h_pT(p, T);

		return h1 + (h2 - h1) * x;
	}
}
