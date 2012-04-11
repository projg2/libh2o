/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "boundaries.h"
#include "consts.h"
#include "region1.h"
#include "region2.h"
#include "region3.h"
#include "region4.h"
#include "region5.h"
#include "xmath.h"

enum h2o_region h2o_region_pT(double p, double T) /* [MPa, K] */
{
	if (T < Tmin || T > Tmax || p <= pmin || p > pmax)
		return H2O_REGION_OUT_OF_RANGE;

	else if (T <= Tb13) /* 1 or 2 */
	{
		if (p > h2o_region4_p_T(T))
			return H2O_REGION1;
		else
			return H2O_REGION2;
	}

	else if (T <= Tb25) /* 3 or 2 */
	{
		if (p >= h2o_b23_p_T(T))
			return H2O_REGION3;
		else
			return H2O_REGION2;
	}

	else /* 5? */
	{
		if (p > pmax5)
			return H2O_REGION_OUT_OF_RANGE;
		else
			return H2O_REGION5;
	}
}

enum h2o_region h2o_region_ph(double p, double h) /* [MPa, kJ/kg] */
{
	if (p < pmin || p > pmax)
		return H2O_REGION_OUT_OF_RANGE;

	/* Left boundary. */
	if (h < h2o_region1_h_pT(p, Tmin))
		return H2O_REGION_OUT_OF_RANGE;

	/* Check the region4 curves. */
	if (p <= psatmax)
	{
		double Tsat = h2o_region4_T_p(p);

		if (h <= h2o_region1_h_pT(p, Tsat))
			return H2O_REGION1;
		else if (h < h2o_region2_h_pT(p, Tsat))
			return H2O_REGION4;
	}
	else /* Then, check the B13 & B23. */
	{
		/* psat3(h) validity range */
		if (h >= 1670.858218 && h <= 2563.592004)
		{
			if (p >= h2o_region3_psat_h(h))
				return H2O_REGION3;
			else
				return H2O_REGION_OUT_OF_RANGE; /* XXX: adjust region 4 */
		}
		else if (h <= h2o_region1_h_pT(p, Tb13))
			return H2O_REGION1;
		else if (h < h2o_region2_h_pT(p, h2o_b23_T_p(p)))
			return H2O_REGION3;
	}

	/* Finally, check B25/right border. */
	if (h <= h2o_region2_h_pT(p, Tb25))
		return H2O_REGION2;
	else if (p <= pmax5 && h <= h2o_region5_s_pT(p, Tmax))
		return H2O_REGION5;
	else
		return H2O_REGION_OUT_OF_RANGE;
}

enum h2o_region h2o_region_ps(double p, double s) /* [MPa, kJ/kgK] */
{
	if (p < pmin || p > pmax || s < 0)
		return H2O_REGION_OUT_OF_RANGE;

	/* First, check the region4 curves. */
	if (p <= psatmax)
	{
		double Tsat = h2o_region4_T_p(p);

		if (s <= h2o_region1_s_pT(p, Tsat))
			return H2O_REGION1;
		else if (s < h2o_region2_s_pT(p, Tsat))
			return H2O_REGION4;
	}
	else /* Then, check the B13 & B23. */
	{
		if (s <= h2o_region1_s_pT(p, Tb13))
			return H2O_REGION1;
		/* psat3(s) validity range */
		else if (s >= 3.778281340 && s <= 5.210887825 && p < h2o_region3_psat_s(s))
			return H2O_REGION_OUT_OF_RANGE; /* XXX: adjust region 4 */
		else if (s < h2o_region2_s_pT(p, h2o_b23_T_p(p)))
			return H2O_REGION3;
	}

	/* Finally, check B25/right border. */
	if (s <= h2o_region2_s_pT(p, Tb25))
		return H2O_REGION2;
	else if (p <= pmax5 && s <= h2o_region5_s_pT(p, Tmax))
		return H2O_REGION5;
	else
		return H2O_REGION_OUT_OF_RANGE;
}

enum h2o_region h2o_region_hs(double h, double s)
{
	/* XXX: check range better */

	if (s < -1.545495919E-4)
		return H2O_REGION_OUT_OF_RANGE;
	else if (s <= 3.778281340) /* B14 */
	{
		if (h < h2o_b14_h_s(s))
			return H2O_REGION4;
		else if (s <= 3.397782955 && h > h2o_b13_h_s(s))
			return H2O_REGION3;
		else
			return H2O_REGION1;
	}
	else if (s <= 4.41202148223476) /* B3a4 */
	{
		if (h < h2o_b3a4_h_s(s))
			return H2O_REGION4;
		else
			return H2O_REGION3;
	}
	else if (s <= 5.85) /* B2c3b4 */
	{
		if (h < h2o_b2c3b4_h_s(s))
			return H2O_REGION4;
		else if (s <= 5.048096828 && s <= 5.260578707)
			return H2O_REGION3;
		else if (s >= 5.260578707 || h >= 2812.942061E3)
			return H2O_REGION2;
		else
		{
			double p = h2o_region2c_p_hs(h, s);

			if (p >= h2o_b23_p_T(h2o_b23_T_hs(h, s)))
				return H2O_REGION3;
			else
				return H2O_REGION2;
		}
	}
	else if (s <= 9.155759395) /* B2ab4 */
	{
		if (h < h2o_b2ab4_h_s(s))
			return H2O_REGION4;
		else
			return H2O_REGION2;
	}
	else /* XXX: min boundary */
		return H2O_REGION2;
}

enum h2o_region h2o_region_Tx(double T, double x) /* [K, 0..1] */
{
	if (x < 0 || x > 1 || T < Tmin || T > Tb13)
		return H2O_REGION_OUT_OF_RANGE;

	return H2O_REGION4;
}

enum h2o_region h2o_region_px(double p, double x) /* [MPa, 0..1] */
{
	if (x < 0 || x > 1 || p < 611.213E-6 || p > psatmax)
		return H2O_REGION_OUT_OF_RANGE;

	return H2O_REGION4;
}

enum h2o_region h2o_region_rhoT(double rho, double T) /* [kg/m³, K] */
{
	/* XXX: support other regions */
	if (T >= Tb13)
	{
		double v = 1/rho;

		if (v <= h2o_region2_v_pT(h2o_b23_p_T(T), T))
			return H2O_REGION3;
		/* XXX: upper limit? */
	}

	return H2O_REGION_OUT_OF_RANGE;
}
