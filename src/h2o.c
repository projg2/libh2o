/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "consts.h"
#include "h2o.h"
#include "region1.h"
#include "region2.h"
#include "region3.h"
#include "region4.h"
#include "region5.h"

#include <assert.h>

typedef double (*twoarg_func_t)(double, double);

h2o_t h2o_new_pT(double p, double T)
{
	h2o_t ret;
	enum h2o_region region = h2o_region_pT(p, T);

	switch (region)
	{
		case H2O_REGION_OUT_OF_RANGE:
			break;
		case H2O_REGION3: /* -> (rho,T) */
			ret._arg1 = 1 / h2o_region3_v_pT(p, T);
			ret._arg2 = T;
			break;
		case H2O_REGION1:
		case H2O_REGION2:
		case H2O_REGION5:
			ret._arg1 = p;
			ret._arg2 = T;
			break;
		default:
			assert(not_reached);
	}
	ret.region = region;

	return ret;
}

h2o_t h2o_new_Tx(double T, double x)
{
	h2o_t ret;
	enum h2o_region region = h2o_region_Tx(T, x);

	switch (region)
	{
		case H2O_REGION_OUT_OF_RANGE:
			break;
		case H2O_REGION4:
			ret._arg1 = T;
			ret._arg2 = x;
			break;
		default:
			assert(not_reached);
	}
	ret.region = region;

	return ret;
}

h2o_t h2o_new_px(double p, double x)
{
	h2o_t ret;
	enum h2o_region region = h2o_region_px(p, x);

	switch (region)
	{
		case H2O_REGION_OUT_OF_RANGE:
			break;
		case H2O_REGION4:
			ret._arg1 = h2o_region4_T_p(p);
			ret._arg2 = x;
			break;
		default:
			assert(not_reached);
	}
	ret.region = region;

	return ret;
}

h2o_t h2o_new_ph(double p, double h)
{
	h2o_t ret;
	enum h2o_region region = h2o_region_ph(p, h);

	switch (region)
	{
		case H2O_REGION4: /* (T, x) */
		{
			double T = h2o_region4_T_p(p);

			ret._arg1 = T;
			ret._arg2 = h2o_region4_x_Th(T, h);
			break;
		}
		case H2O_REGION5:
			region = H2O_REGION_OUT_OF_RANGE;
			break;
		case H2O_REGION_OUT_OF_RANGE:
			break;

		default:
		{
			twoarg_func_t T_getter;

			switch (region)
			{
				case H2O_REGION1:
					T_getter = &h2o_region1_T_ph;
					break;
				case H2O_REGION2:
					T_getter = &h2o_region2_T_ph;
					break;
				case H2O_REGION3:
					T_getter = &h2o_region3_T_ph;
					break;
				default:
					assert(not_reached);
			}

			if (region == H2O_REGION3)
				ret._arg1 = 1 / h2o_region3_v_ph(p, h);
			else
				ret._arg1 = p;

			ret._arg2 = T_getter(p, h);
		}
	}
	ret.region = region;

	return ret;
}

h2o_t h2o_new_ps(double p, double s)
{
	h2o_t ret;
	enum h2o_region region = h2o_region_ps(p, s);

	switch (region)
	{
		case H2O_REGION4: /* (T, x) */
		{
			double T = h2o_region4_T_p(p);

			ret._arg1 = T;
			ret._arg2 = h2o_region4_x_Ts(T, s);
			break;
		}
		case H2O_REGION5:
			region = H2O_REGION_OUT_OF_RANGE;
			break;
		case H2O_REGION_OUT_OF_RANGE:
			break;

		default:
		{
			twoarg_func_t T_getter;

			switch (region)
			{
				case H2O_REGION1:
					T_getter = &h2o_region1_T_ps;
					break;
				case H2O_REGION2:
					T_getter = &h2o_region2_T_ps;
					break;
				case H2O_REGION3:
					T_getter = &h2o_region3_T_ps;
					break;
				default:
					assert(not_reached);
			}

			if (region == H2O_REGION3)
				ret._arg1 = 1 / h2o_region3_v_ps(p, s);
			else
				ret._arg1 = p;

			ret._arg2 = T_getter(p, s);
		}
	}
	ret.region = region;

	return ret;
}

h2o_t h2o_new_hs(double h, double s)
{
	h2o_t ret;
	enum h2o_region region = h2o_region_hs(h, s);

	switch (region)
	{
		case H2O_REGION5:
			ret.region = H2O_REGION_OUT_OF_RANGE;
			break;
		case H2O_REGION_OUT_OF_RANGE:
			ret.region = region;
			break;
		default:
		{
			twoarg_func_t getter;
			double arg1;

			switch (region)
			{
				case H2O_REGION1:
					getter = &h2o_region1_p_hs;
					break;
				case H2O_REGION2:
					getter = &h2o_region2_p_hs;
					break;
				case H2O_REGION3:
					getter = &h2o_region3_p_hs;
					break;
				case H2O_REGION4:
					getter = &h2o_region4_T_hs;
					break;
				default:
					assert(not_reached);
			}

			arg1 = getter(h, s);

			if (region == H2O_REGION4)
				ret = h2o_new_Tx(arg1, h2o_region4_x_Th(arg1, h));
			else
				ret = h2o_new_ps(arg1, s);
		}
	}

	return ret;
}

h2o_t h2o_new_rhoT(double rho, double T)
{
	h2o_t ret;
	enum h2o_region region = h2o_region_rhoT(rho, T);

	switch (region)
	{
		case H2O_REGION_OUT_OF_RANGE:
			break;
		case H2O_REGION3:
			ret._arg1 = rho;
			ret._arg2 = T;
			break;
		default:
			assert(not_reached);
	}
	ret.region = region;

	return ret;
}

int h2o_is_valid(const h2o_t state)
{
	return state.region != H2O_REGION_OUT_OF_RANGE;
}

enum h2o_region h2o_get_region(const h2o_t state)
{
	return state.region;
}

double h2o_get_p(const h2o_t state)
{
	double ret;

	switch (state.region)
	{
		case H2O_REGION1:
		case H2O_REGION2:
		case H2O_REGION5:
			ret = state._arg1;
			break;
		case H2O_REGION3:
			ret = h2o_region3_p_rhoT(state._arg1, state._arg2);
			break;
		case H2O_REGION4:
			ret = h2o_region4_p_T(state._arg1);
			break;
		default:
			assert(not_reached);
	}

	return ret;
}

double h2o_get_T(const h2o_t state)
{
	double ret;

	switch (state.region)
	{
		case H2O_REGION1:
		case H2O_REGION2:
		case H2O_REGION3:
		case H2O_REGION5:
			ret = state._arg2;
			break;
		case H2O_REGION4:
			ret = state._arg1;
			break;
		default:
			assert(not_reached);
	}

	return ret;
}

double h2o_get_x(const h2o_t state)
{
	double ret;

	switch (state.region)
	{
		case H2O_REGION1: /* water */
			ret = 0;
			break;
		case H2O_REGION2:
		case H2O_REGION5: /* dry steam */
			ret = 1;
			break;
		case H2O_REGION4: /* wet steam */
			ret = state._arg2;
			break;
		default:
			assert(not_reached);
	}

	return ret;
}

double h2o_get_rho(const h2o_t state)
{
	double ret;

	switch (state.region)
	{
		case H2O_REGION3:
			ret = state._arg1;
			break;
		default:
			ret = 1 / h2o_get_v(state);
	}

	return ret;
}

static double region3_v_rhoT(double rho, double T)
{
	return 1 / rho;
}

double h2o_get_v(const h2o_t state)
{
	twoarg_func_t func;

	switch (state.region)
	{
		case H2O_REGION1:
			func = &h2o_region1_v_pT;
			break;
		case H2O_REGION2:
			func = &h2o_region2_v_pT;
			break;
		case H2O_REGION3:
			func = &region3_v_rhoT;
			break;
		case H2O_REGION4:
			func = &h2o_region4_v_Tx;
			break;
		case H2O_REGION5:
			func = &h2o_region5_v_pT;
			break;
		default:
			assert(not_reached);
	}

	return func(state._arg1, state._arg2);
}

double h2o_get_u(const h2o_t state)
{
	twoarg_func_t func;

	switch (state.region)
	{
		case H2O_REGION1:
			func = &h2o_region1_u_pT;
			break;
		case H2O_REGION2:
			func = &h2o_region2_u_pT;
			break;
		case H2O_REGION3:
			func = &h2o_region3_u_rhoT;
			break;
		case H2O_REGION4:
			func = &h2o_region4_u_Tx;
			break;
		case H2O_REGION5:
			func = &h2o_region5_u_pT;
			break;
		default:
			assert(not_reached);
	}

	return func(state._arg1, state._arg2);
}

double h2o_get_h(const h2o_t state)
{
	twoarg_func_t func;

	switch (state.region)
	{
		case H2O_REGION1:
			func = &h2o_region1_h_pT;
			break;
		case H2O_REGION2:
			func = &h2o_region2_h_pT;
			break;
		case H2O_REGION3:
			func = &h2o_region3_h_rhoT;
			break;
		case H2O_REGION4:
			func = &h2o_region4_h_Tx;
			break;
		case H2O_REGION5:
			func = &h2o_region5_h_pT;
			break;
		default:
			assert(not_reached);
	}

	return func(state._arg1, state._arg2);
}

double h2o_get_s(const h2o_t state)
{
	twoarg_func_t func;

	switch (state.region)
	{
		case H2O_REGION1:
			func = &h2o_region1_s_pT;
			break;
		case H2O_REGION2:
			func = &h2o_region2_s_pT;
			break;
		case H2O_REGION3:
			func = &h2o_region3_s_rhoT;
			break;
		case H2O_REGION4:
			func = &h2o_region4_s_Tx;
			break;
		case H2O_REGION5:
			func = &h2o_region5_s_pT;
			break;
		default:
			assert(not_reached);
	}

	return func(state._arg1, state._arg2);
}

double h2o_get_cp(const h2o_t state)
{
	twoarg_func_t func;

	switch (state.region)
	{
		case H2O_REGION1:
			func = &h2o_region1_cp_pT;
			break;
		case H2O_REGION2:
			func = &h2o_region2_cp_pT;
			break;
		case H2O_REGION3:
			func = &h2o_region3_cp_rhoT;
			break;
		case H2O_REGION4:
			func = &h2o_region4_cp_Tx;
			break;
		case H2O_REGION5:
			func = &h2o_region5_cp_pT;
			break;
		default:
			assert(not_reached);
	}

	return func(state._arg1, state._arg2);
}

double h2o_get_cv(const h2o_t state)
{
	twoarg_func_t func;

	switch (state.region)
	{
		case H2O_REGION1:
			func = &h2o_region1_cv_pT;
			break;
		case H2O_REGION2:
			func = &h2o_region2_cv_pT;
			break;
		case H2O_REGION3:
			func = &h2o_region3_cv_rhoT;
			break;
		case H2O_REGION4:
			func = &h2o_region4_cv_Tx;
			break;
		case H2O_REGION5:
			func = &h2o_region5_cv_pT;
			break;
		default:
			assert(not_reached);
	}

	return func(state._arg1, state._arg2);
}

double h2o_get_w(const h2o_t state)
{
	twoarg_func_t func;

	switch (state.region)
	{
		case H2O_REGION1:
			func = &h2o_region1_w_pT;
			break;
		case H2O_REGION2:
			func = &h2o_region2_w_pT;
			break;
		case H2O_REGION3:
			func = &h2o_region3_w_rhoT;
			break;
		case H2O_REGION4:
			func = &h2o_region4_w_Tx;
			break;
		case H2O_REGION5:
			func = &h2o_region5_w_pT;
			break;
		default:
			assert(not_reached);
	}

	return func(state._arg1, state._arg2);
}
