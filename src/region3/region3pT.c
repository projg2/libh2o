/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <stdlib.h>

#include "consts.h"
#include "region3.h"
#include "region4.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Specific Volume
 * as a Function of Pressure and Temperature v(p,T)
 * for Region 3 of the IAPWS Industrial Formulation 1997 for the
 * Thermodynamic Properties of Water and Steam */

enum h2o_region3_subregion_pT h2o_region3_subregion_pT(double p, double T)
{
	if (p > 40)
	{
		if (T <= h2o_region3ab_T_p(p))
			return H2O_REGION3A_PT;
		else
			return H2O_REGION3B_PT;
	}
	else if (p <= p3cd)
	{
		if (T <= h2o_region4_T_p(p))
			return H2O_REGION3C_PT;
		else
			return H2O_REGION3T_PT;
	}
	else if (T <= h2o_region3cd_T_p(p))
		return H2O_REGION3C_PT;
	else if (p > 25)
	{
		if (T <= h2o_region3ab_T_p(p))
			return H2O_REGION3D_PT;
		else if (T <= h2o_region3ef_T_p(p))
			return H2O_REGION3E_PT;
		else
			return H2O_REGION3F_PT;
	}
	else if (p <= 20.5)
	{
		if (T <= h2o_region4_T_p(p))
			return H2O_REGION3S_PT;
		else
			return H2O_REGION3T_PT;
	}
	else if (T > h2o_region3jk_T_p(p))
		return H2O_REGION3K_PT;
	else if (p > 23)
	{
		if (T <= h2o_region3gh_T_p(p))
		{
			if (p > 23.5)
				return H2O_REGION3G_PT;
			else
				return H2O_REGION3L_PT;
		}
		else if (T <= h2o_region3ef_T_p(p))
			return H2O_REGION3H_PT;
		else if (T <= h2o_region3ij_T_p(p))
			return H2O_REGION3I_PT;
		else /*if (T <= h2o_region3jk_T_p(p))*/
			return H2O_REGION3J_PT;
	}
	else if (p > 22.5)
	{
		if (T <= h2o_region3gh_T_p(p))
			return H2O_REGION3L_PT;
		else if (T <= h2o_region3mn_T_p(p))
			return H2O_REGION3M_PT;
		else if (T <= h2o_region3ef_T_p(p))
			return H2O_REGION3N_PT;
		else if (T <= h2o_region3op_T_p(p))
			return H2O_REGION3O_PT;
		else if (T <= h2o_region3ij_T_p(p))
			return H2O_REGION3P_PT;
		else /*if (T <= h2o_region3jk_T_p(p))*/
			return H2O_REGION3J_PT;
	}
	else if (p > psat3rsmax) /* @ 643.15K */
	{
		if (T <= h2o_region3qu_T_p(p))
			return H2O_REGION3Q_PT;
		else if (T <= h2o_region3rx_T_p(p))
		{
			/* backwards equations */
			if (p > 22.11)
			{
				if (T <= h2o_region3uv_T_p(p))
					return H2O_REGION3U_PT;
				else if (T <= h2o_region3ef_T_p(p))
					return H2O_REGION3V_PT;
				else if (T <= h2o_region3wx_T_p(p))
					return H2O_REGION3W_PT;
				else
					return H2O_REGION3X_PT;
			}
			else if (p > pcrit)
			{
				if (T <= h2o_region3uv_T_p(p))
					return H2O_REGION3U_PT;
				else if (T <= h2o_region3ef_T_p(p))
					return H2O_REGION3Y_PT;
				else if (T <= h2o_region3wx_T_p(p))
					return H2O_REGION3Z_PT;
				else
					return H2O_REGION3X_PT;
			}
			else
			{
				if (T <= h2o_region4_T_p(p))
				{
					if (p > p3ymin && T > h2o_region3uv_T_p(p))
						return H2O_REGION3Y_PT;
					else
						return H2O_REGION3U_PT;
				}
				else
				{
					if (p > p3zmin && T <= h2o_region3wx_T_p(p))
						return H2O_REGION3Z_PT;
					else
						return H2O_REGION3X_PT;
				}
			}
		}
		else /*if (T <= h2o_region3jk_T_p(p))*/
			return H2O_REGION3R_PT;
	}
	else /*if (p > 20.5)*/
	{
		if (T <= h2o_region4_T_p(p))
			return H2O_REGION3S_PT;
		else /*if (T <= h2o_region3jk_T_p(p))*/
			return H2O_REGION3R_PT;
	}
}

typedef double (*twoarg_func_t)(double, double);

/* this function intends to make compiler happy. */
static double impossible_happened(double a, double b)
{
	abort();

	if (a || b)
	{
	}
}

double h2o_region3_v_pT(double p, double T)
{
	twoarg_func_t v_getter;

	switch (h2o_region3_subregion_pT(p, T))
	{
		case H2O_REGION3A_PT:
			v_getter = &h2o_region3a_v_pT;
			break;
		case H2O_REGION3B_PT:
			v_getter = &h2o_region3b_v_pT;
			break;
		case H2O_REGION3C_PT:
			v_getter = &h2o_region3c_v_pT;
			break;
		case H2O_REGION3D_PT:
			v_getter = &h2o_region3d_v_pT;
			break;
		case H2O_REGION3E_PT:
			v_getter = &h2o_region3e_v_pT;
			break;
		case H2O_REGION3F_PT:
			v_getter = &h2o_region3f_v_pT;
			break;
		case H2O_REGION3G_PT:
			v_getter = &h2o_region3g_v_pT;
			break;
		case H2O_REGION3H_PT:
			v_getter = &h2o_region3h_v_pT;
			break;
		case H2O_REGION3I_PT:
			v_getter = &h2o_region3i_v_pT;
			break;
		case H2O_REGION3J_PT:
			v_getter = &h2o_region3j_v_pT;
			break;
		case H2O_REGION3K_PT:
			v_getter = &h2o_region3k_v_pT;
			break;
		case H2O_REGION3L_PT:
			v_getter = &h2o_region3l_v_pT;
			break;
		case H2O_REGION3M_PT:
			v_getter = &h2o_region3m_v_pT;
			break;
		case H2O_REGION3N_PT:
			v_getter = &h2o_region3n_v_pT;
			break;
		case H2O_REGION3O_PT:
			v_getter = &h2o_region3o_v_pT;
			break;
		case H2O_REGION3P_PT:
			v_getter = &h2o_region3p_v_pT;
			break;
		case H2O_REGION3Q_PT:
			v_getter = &h2o_region3q_v_pT;
			break;
		case H2O_REGION3R_PT:
			v_getter = &h2o_region3r_v_pT;
			break;
		case H2O_REGION3S_PT:
			v_getter = &h2o_region3s_v_pT;
			break;
		case H2O_REGION3T_PT:
			v_getter = &h2o_region3t_v_pT;
			break;
		case H2O_REGION3U_PT:
			v_getter = &h2o_region3u_v_pT;
			break;
		case H2O_REGION3V_PT:
			v_getter = &h2o_region3v_v_pT;
			break;
		case H2O_REGION3W_PT:
			v_getter = &h2o_region3w_v_pT;
			break;
		case H2O_REGION3X_PT:
			v_getter = &h2o_region3x_v_pT;
			break;
		case H2O_REGION3Y_PT:
			v_getter = &h2o_region3y_v_pT;
			break;
		case H2O_REGION3Z_PT:
			v_getter = &h2o_region3z_v_pT;
			break;
		default:
			v_getter = &impossible_happened;
	}

	return v_getter(p, T);
}
