/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region2.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Pressure as a Function
 * of Enthalpy and Entropy p(h,s) to the IAPWS Industrial Formulation 1997
 * for the Thermodynamic Properties of Water and Steam;
 * s. 6: Backward Equation p(h,s) for Region 1 */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	+0.801496989929495E-01, -0.543862807146111E+00,
	+0.337455597421283E+00, +0.890555451157450E+01,
	+0.313840736431485E+03, +0.797367065977789E+00,
	-0.121616973556240E+01, +0.872803386937477E+01,

	-0.169769781757602E+02, -0.186552827328416E+03,
	+0.951159274344237E+05, -0.189168510120494E+02,
	-0.433407037194840E+04, +0.543212633012715E+09,
	+0.144793408386013E+00, +0.128024559637516E+03,

	-0.672309534071268E+05, +0.336972380095287E+08,
	-0.586634196762720E+03, -0.221403224769889E+11,
	+0.171606668708389E+04, -0.570817595806302E+09,
	-0.312109693178482E+04, -0.207841384633010E+07,

	+0.305605946157786E+13, +0.322157004314333E+04,
	+0.326810259797295E+12, -0.144104158934487E+04,
	+0.410694867802691E+03, +0.109077066873024E+12,
	-0.247964654258893E+14, +0.188801906865134E+10,

	-0.123651009018773E+15
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 14
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0,
	1, 1, 1, 1, 1, 1,
	2, 2, 2, 3, 3, 3, 3,
	4, 4, 5, 5, 6, 6, 6,
	7, 7, 8, 8, 8, 8, 9, 10
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18
};

static const int J[] = {
	0,

	0, 1, 2, 4, 8,
	0, 1, 2, 3, 5, 10,
	1, 6, 13, 0, 1, 7, 10,
	1, 12, 1, 10, 1, 8, 13,
	1, 12, 1, 3, 11, 13, 9, 12
};

static const double pstar = 100; /* [MPa] */
static const double hstar = 4100; /* [kJ/kg] */
static const double sstar = 7.9; /* [kJ/kgK] */

double h2o_region2b_p_hs(double h, double s)
{
	double eta = h / hstar;
	double sigma = s / sstar;

	double sum = poly_value(eta - 0.6, sigma - 1.01,
			I, Ipows, 0, 11, 0,
			J, Jpows, 0, 14, 0,
			n, 33);

	return pow4(sum) * pstar;
}
