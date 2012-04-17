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

#include <assert.h>
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

static inline double poly_value(double x,
		int min, int max, int deriv,
		const double n[])
{
	double sum;

	const int maxn = 1 - min + max;

	int i;

	sum = n[maxn];

	if (deriv == 1)
		sum *= maxn + min - 1;

	for (i = maxn - 1; i >= 1; --i)
	{
		double coeff = n[i];

		if (deriv == 1)
			coeff *= i + min - 1;

		sum *= x;
		sum += coeff;
	}

	if (min - deriv != 0)
		sum *= pow(x, min - deriv);

	return sum;
}

static inline double twoarg_poly_value(double x1, double x2,
		const int I[], const double Ipows[], int Ipowzero,
		int Ipowlen, int x1der,
		const int J[], const double Jpows[], int Jpowzero,
		int Jpowlen, int x2der,
		const double n[], int nlen)
{
	double sum = 0;

	int i;

	double x1powers[20], x2powers[34];

	assert(Ipowlen <= 20);
	assert(Jpowlen <= 34);

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
