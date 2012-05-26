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

static inline double pow2(double arg);
static inline double pow4(double arg);
static inline double quadr_value(double a, double b, double c, double x);

double poly_value(double x,
		int min, int max, int deriv,
		const double n[]);
double twoarg_poly_value(double x1, double x2,
		const int I[], const double Ipows[], int Ipowzero,
		int Ipowlen, int x1der,
		const int J[], const double Jpows[], int Jpowzero,
		int Jpowlen, int x2der,
		const double n[], int nlen);

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

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_XMATH_H*/
