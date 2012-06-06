/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "consts.h"
#include "region2.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 6: Equations for Region 2 */

/* ideal-gas part coefficients */
static const double no[] = {
	+0.00000000000000E+0,

	-0.56087911283020E-2, +0.71452738081455E-1,
	-0.40710498223928E+0, +0.14240819171444E+1,
	-0.43839511319450E+1, -0.96937268393049E+1,
	+0.10087275970006E+2, -0.28408632460772E+0,
	+0.21268463753307E-1
};

/* resident part coefficients */
static const double n[] = {
	+0.00000000000000E+00,

	-0.73362260186506E-2, -0.88223831943146E-1,
	-0.72334555213245E-1, -0.40813178534455E-2,
	+0.20097803380207E-2, -0.53045921898642E-1,
	-0.76190409086970E-2, -0.63498037657313E-2,

	-0.86043093028588E-1, +0.75321581522770E-2,
	-0.79238375446139E-2, -0.22888160778447E-3,
	-0.26456501482810E-2
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5
};

static const int I[] = {
	0,

	1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5
};

static const double Jpows[] = {
	0, 1, 2, 4, 5, 7, 9, 10, 11, 16
};

static const int J[] = {
	0,

	0, 2, 4, 8, 1, 5, 9, 3, 9, 5, 7, 6, 7
};

static const double Tstar = 540; /* [K] */

static double h2o_region2_meta_gammao_pitau(double pi, double tau, int pider, int tauder)
{
	if (!pider)
	{
		double sum = poly_value(tau, -5, 3, tauder, no);

		if (!tauder)
			sum += log(pi);

		return sum;
	}
	else
		return 1/pi;
}

static double h2o_region2_meta_gammar_pitau(double pi, double tau, int pider, int tauder)
{
	return twoarg_poly_value(pi, tau - 0.5,
			I, Ipows, 0, 6, pider,
			J, Jpows, 0, 10, tauder,
			n, 13);
}

static double h2o_region2_meta_gamma_pitau(double pi, double tau, int pider, int tauder)
{
	double sum;

	sum = h2o_region2_meta_gammao_pitau(pi, tau, pider, tauder);
	sum += h2o_region2_meta_gammar_pitau(pi, tau, pider, tauder);

	return sum;
}

double h2o_region2_meta_v_pT(double p, double T)
{
	double pi = p;
	double tau = Tstar / T;

	double gammapi = h2o_region2_meta_gamma_pitau(pi, tau, 1, 0);

	return pi * gammapi * R * T / p * 1E-3;
}

double h2o_region2_meta_u_pT(double p, double T)
{
	double pi = p;
	double tau = Tstar / T;

	double gammatau = h2o_region2_meta_gamma_pitau(pi, tau, 0, 1);
	double gammapi = h2o_region2_meta_gamma_pitau(pi, tau, 1, 0);

	return (tau * gammatau - pi * gammapi) * R * T;
}

double h2o_region2_meta_s_pT(double p, double T)
{
	double pi = p;
	double tau = Tstar / T;

	double gammatau = h2o_region2_meta_gamma_pitau(pi, tau, 0, 1);
	double gamma = h2o_region2_meta_gamma_pitau(pi, tau, 0, 0);

	return (tau * gammatau - gamma) * R;
}

double h2o_region2_meta_h_pT(double p, double T)
{
	double pi = p;
	double tau = Tstar / T;

	double gammatau = h2o_region2_meta_gamma_pitau(pi, tau, 0, 1);

	return tau * gammatau * R * T;
}
