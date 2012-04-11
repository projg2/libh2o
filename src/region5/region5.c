/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif


#include "consts.h"
#include "region5.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 9: Equations for Region 5 */

/* ideal-gas part coefficients */
static const double no_store[] = {
	-0.24805148933466E-1,
	+0.36901534980333E+0,
	-0.31161318213925E+1,
	-0.13179983674201E+2,
	+0.68540841634434E+1,
	-0.32961626538917E+0
};

/* shift for negative indices */
static const double *no = &no_store[3];

/* resident part coefficients */
static const double n[] = {
	+0.00000000000000E+00,

	+0.15736404855259E-2, +0.90153761673944E-3,
	-0.50270077677648E-2, +0.22440037409485E-5,
	-0.41163275453471E-5, +0.37919454822955E-7
};

static const double Ipows[] = {
	0, 1, 2, 3
};

static const int I[] = {
	0,

	1, 1, 1, 2, 2, 3
};

static const double Jpows[] = {
	0, 1, 2, 3, 7, 9
};

static const int J[] = {
	0,

	1, 2, 3, 3, 5, 4
};

static const double Tstar = 1000; /* [K] */

static inline double h2o_region5_gammao_pitau(double pi, double tau, int pider, int tauder)
{
	if (!pider)
	{
		int i;
		double sum = 0;

		double taupowers_store[3 + 3];
			/* shift it for negative indices */
		double* taupowers = &taupowers_store[3];

		fill_powers_incr(taupowers, 3, tau, tauder);
		fill_powers_decr(taupowers, -3, tau, tauder);

		if (!tauder)
			sum = log(pi);

		for (i = -3; i <= 2; ++i)
		{
			double taupow = taupowers[i];

			double memb = no[i] * taupow;
			if (tauder == 1)
				memb *= i;

			sum += memb;
		}

		return sum;
	}
	else
		return 1/pi;
}

static inline double h2o_region5_gammar_pitau(double pi, double tau, int pider, int tauder)
{

	return poly_value(pi, tau,
			I, Ipows, 0, 4, pider,
			J, Jpows, 0, 6, tauder,
			n, 6);
}

static inline double h2o_region5_gamma_pT(double p, double T, int pider, int tauder)
	/* p [MPa], T [K], pider, tauder: 0/1 */
{
	double tau = Tstar / T;

	double sum;

	/* ideal-gas part */
	sum = h2o_region5_gammao_pitau(p, tau, pider, tauder);
	sum += h2o_region5_gammar_pitau(p, tau, pider, tauder);

	if (pider == 1)
		sum *= p;
	if (tauder == 1)
		sum *= tau;

	return sum;
}

double h2o_region5_v_pT(double p, double T) /* [MPa, K] -> [m³/kg] */
{
	double gammapi = h2o_region5_gamma_pT(p, T, 1, 0);

	return gammapi * R * T / p * 1E-3;
}

double h2o_region5_u_pT(double p, double T) /* [MPa, K] -> [kJ/kg] */
{
	double gammatau = h2o_region5_gamma_pT(p, T, 0, 1);
	double gammapi = h2o_region5_gamma_pT(p, T, 1, 0);

	return (gammatau - gammapi) * R * T;
}

double h2o_region5_s_pT(double p, double T) /* [MPa, K] -> [kJ/kgK] */
{
	double gammatau = h2o_region5_gamma_pT(p, T, 0, 1);
	double gamma = h2o_region5_gamma_pT(p, T, 0, 0);

	return (gammatau - gamma) * R;
}

double h2o_region5_h_pT(double p, double T) /* [MPa, K] -> [kJ/kg] */
{
	double gammatau = h2o_region5_gamma_pT(p, T, 0, 1);

	return gammatau * R * T;
}
