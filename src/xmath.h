/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_XMATH_H
#define _H2O_XMATH_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

#include <math.h>

/* using pow(arg, N) is basically equivalent
 * but it will be optimized out only with -ffast-math;
 * this way, optimization always takes place. */

static inline double pow2(double arg)
{
	return arg * arg;
}

static inline double pow4(double arg)
{
	return pow2(pow2(arg));
}

static inline double quadr_value(double a, double b, double c, double x)
{
	return (a * x + b) * x + c;
}

static inline double deg3_value(double a, double b, double c, double d, double x)
{
	return quadr_value(a, b, c, x) * x + d;
}

static inline void fill_powers(double* powers, const double* exponents,
		int zeropos, int count, double expr, int deriv)
{
	int i;
	double tmp = expr;

	powers[zeropos + deriv] = 1;
	powers[zeropos + deriv + 1] = tmp;

	for (i = zeropos + deriv + 2; i < count; ++i)
	{
		if (exponents[i] - 1 == exponents[i - 1])
			tmp *= expr;
		else
			tmp = pow(expr, exponents[i] - deriv);
		powers[i] = tmp;
	}

	tmp = 1.0;
	for (i = zeropos + deriv - 1; i >= 0; --i)
	{
		if (exponents[i] + 1 == exponents[i + 1])
			tmp /= expr;
		else
			tmp = pow(expr, exponents[i] - deriv);
		powers[i] = tmp;
	}
}

static inline void fill_powers_incr(double* powers, int count, double expr, int deriv)
{
	int i;
	double tmp = expr;

	powers[deriv] = 1;
	powers[deriv + 1] = tmp;

	for (i = deriv + 2; i < count; ++i)
	{
		tmp *= expr;
		powers[i] = tmp;
	}

	if (deriv)
		powers[0] = 1/expr;
}

static inline void fill_powers_decr(double* powers, int min, double expr, int deriv)
{
	int i;
	double tmp = 1.0;

	if (deriv)
		tmp /= expr;

	for (i = -1; i >= min; --i)
	{
		tmp /= expr;
		powers[i] = tmp;
	}
}

static inline double poly_value(double x1, double x2,
		const int I[], const double Ipows[], int Ipowzero,
		int Ipowlen, int x1der,
		const int J[], const double Jpows[], int Jpowzero,
		int Jpowlen, int x2der,
		const double n[], int nlen)
{
	double sum = 0;

	int i;

	double x1powers[Ipowlen], x2powers[Jpowlen];

	fill_powers(x1powers, Ipows, Ipowzero, Ipowlen, x1, x1der);
	fill_powers(x2powers, Jpows, Jpowzero, Jpowlen, x2, x2der);

	for (i = 1; i <= nlen; ++i)
	{
		double x1pow = x1powers[I[i]];
		double x2pow = x2powers[J[i]];

		double memb = n[i] * x1pow * x2pow;
		if (x1der == 1)
			memb *= Ipows[I[i]];
		if (x2der == 1)
			memb *= Jpows[J[i]];

		sum += memb;
	}

	return sum;
}

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_XMATH_H*/
