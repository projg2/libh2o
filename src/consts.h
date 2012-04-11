/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_CONSTS_H
#define _H2O_CONSTS_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

/* Based on IF97-Rev, s. 3: Reference Constants */

static const double R = 0.461526; /* [kJ/kgK] */

/* critical point */

static const double Tcrit = 647.096; /* [K] */
static const double pcrit = 22.064; /* [MPa] */
static const double rhocrit = 322; /* [kg/m³] */

/* region boundaries f(p, T) */

static const double Tmin = 273.15;
static const double Tb13 = 623.15;
static const double Tb25 = 1073.15;
static const double Tmax = 2273.15;

static const double pmin = 0;
static const double pmax = 100;
static const double pmax5 = 50;

/* max psat handled by Region 4 */
static const double psatmax = 16.5291642;

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_CONSTS_H*/
