/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif


#include "consts.h"
#include "region3.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 8: Equations for Region 1 */

/* coefficient table; n[0] added for convenience, n[1] is for ln */
static const double n[] = {
	+0.00000000000000E+00,

	+0.10658070028513E+1, -0.15732845290239E+2,
	+0.20944396974307E+2, -0.76867707878716E+1,
	+0.26185947787954E+1, -0.28080781148620E+1,
	+0.12053369696517E+1, -0.84566812812502E-2,

	-0.12654315477714E+1, -0.11524407806681E+1,
	+0.88521043984318E+0, -0.64207765181607E+0,
	+0.38493460186671E+0, -0.85214708824206E+0,
	+0.48972281541877E+1, -0.30502617256965E+1,

	+0.39420536879154E-1, +0.12558408424308E+0,
	-0.27999329698710E+0, +0.13899799569460E+1,
	-0.20189915023570E+1, -0.82147637173963E-2,
	-0.47596035734923E+0, +0.43984074473500E-1,

	-0.44476435428739E+0, +0.90572070719733E+0,
	+0.70522450087967E+0, +0.10770512626332E+0,
	-0.32913623258954E+0, -0.50871062041158E+0,
	-0.22175400873096E-1, +0.94260751665092E-1,

	+0.16436278447961E+0, -0.13503372241348E-1,
	-0.14834345352472E-1, +0.57922953628084E-3,
	+0.32308904703711E-2, +0.80964802996215E-4,
	-0.16557679795037E-3, -0.44923899061815E-4
};

static const int I[] = {
	0, 0,

	0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
	2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
	4, 4, 4, 4, 5, 5, 5, 6, 6, 6,
	7, 8, 9, 9, 10, 10, 11
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 6, 7, 10, 12, 15, 16, 17, 22, 23, 26
};

/* using Jpows[] indexes */
static const int J[] = {
	0, 0,

	0, 1, 2, 6, 7, 8, 13, 2, 5, 9, 11,
	0, 2, 5, 6, 12, 14, 0, 2, 4, 10, 14,
	0, 2, 4, 14, 1, 3, 14, 0, 2, 14,
	2, 14, 2, 14, 0, 1, 14
};

static const double rhostar = rhocrit;
static const double pstar = pcrit;
static const double Tstar = Tcrit;

static inline double h2o_region3_phi_rhoT(double rho, double T, int deltader, int tauder)
	/* rho [kg/m³], T [K], pider, tauder: 0/1 */
{
	double delta = rho / rhostar;
	double tau = Tstar / T;

	double sum;

	int i;

	double deltapowers[12], taupowers[15];

	fill_powers_incr(deltapowers, 12, delta, deltader);
	fill_powers(taupowers, Jpows, 0, 15, tau, tauder);

	if (tauder)
		sum = 0;
	else if (deltader)
		sum = n[1] / delta;
	else
		sum = n[1] * log(delta);

	for (i = 2; i <= 40; ++i)
	{
		double deltapow = deltapowers[I[i]];
		double taupow = taupowers[J[i]];

		double memb = n[i] * deltapow * taupow;
		if (deltader == 1)
			memb *= I[i];
		if (tauder == 1)
			memb *= Jpows[J[i]];

		sum += memb;
	}

	if (deltader == 1)
		sum *= delta;
	if (tauder == 1)
		sum *= tau;

	return sum;
}

double h2o_region3_p_rhoT(double rho, double T) /* [kg/m³, K] -> [MPa] */
{
	double phidelta = h2o_region3_phi_rhoT(rho, T, 1, 0);

	return phidelta * rho * R * T * 1E-3;
}

double h2o_region3_u_rhoT(double rho, double T) /* [kg/m³, K] -> [kJ/kg] */
{
	double phitau = h2o_region3_phi_rhoT(rho, T, 0, 1);

	return phitau * R * T;
}

double h2o_region3_s_rhoT(double rho, double T) /* [kg/m³, K] -> [kJ/kgK] */
{
	double phitau = h2o_region3_phi_rhoT(rho, T, 0, 1);
	double phi = h2o_region3_phi_rhoT(rho, T, 0, 0);

	return (phitau - phi) * R;
}

double h2o_region3_h_rhoT(double rho, double T) /* [kg/m³, K] -> [kJ/kg] */
{
	double phidelta = h2o_region3_phi_rhoT(rho, T, 1, 0);
	double phitau = h2o_region3_phi_rhoT(rho, T, 0, 1);

	return (phitau + phidelta) * R * T;
}
