/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Revised Supplementary Release on Backward Equations for the Functions
 * T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial
 * Formulation 1997 for the Thermodynamic Properties of Water and Steam
 * s. 4.3: Boundary Equations psat(h) and psat(s) */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E00,

	+0.600073641753024E00, -0.936203654849857E01,
	+0.246590798594147E02, -0.107014222858224E03,
	-0.915821315805768E14, -0.862332011700662E04,
	-0.235837344740032E02, +0.252304969384128E18,

	-0.389718771997719E19, -0.333775713645296E23,
	+0.356499469636328E11, -0.148547544720641E27,
	+0.330611514838798E19, +0.813641294467829E38
};

static const double Ipows[] = {
	0, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36
};

static const int I[] = {
	0,

	0, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
};

static const double Jpows[] = {
	0, 1, 3, 4, 8, 16, 18, 24, 36
};

static const int J[] = {
	0,

	0, 1, 2, 3, 8, 2, 0, 7, 5, 5, 2, 6, 4, 7
};

static const double pstar = 22; /* [MPa] */
static const double hstar = 2600; /* [kJ/kg] */
	
double h2o_region3_psat_h(double h)
{
	double eta = h / hstar;
	double etaexprI = eta - 1.02;
	double etaexprJ = eta - 0.608;

	double sum = 0;

	int i;

	double etapowersI[11], etapowersJ[9];

	fill_powers(etapowersI, Ipows, 0, 11, etaexprI, 0);
	fill_powers(etapowersJ, Jpows, 0, 9, etaexprJ, 0);

	for (i = 1; i <= 14; ++i)
	{
		double etapowI = etapowersI[I[i]];
		double etapowJ = etapowersJ[J[i]];

		sum += n[i] * etapowI * etapowJ;
	}

	return sum * pstar;
}
