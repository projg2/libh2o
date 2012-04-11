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
static const double scrit = 4.41202148223476; /* [kJ/kgK] */

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
static const double psatmin = 611.213E-6;

/* Region 3 psat equations validity range */
static const double psat3_hmin = 1670.858218;
static const double psat3_hmax = 2563.592004;
static const double psat3_smin = 3.778281340;
static const double psat3_smax = 5.210887825;

/* region boundaries f(h, s) */

static const double smin = -1.545495919E-4; /* XXX? */
static const double s1max = 3.778281340;
static const double s3min = 3.397782955;
static const double s2cmax = 5.85;
static const double s4max = 9.155759395;

static const double sb23min = 5.048096828;
static const double sb23max = 5.260578707;
static const double hb23min = 2563.592004;
static const double hb23max = 2812.942061;

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_CONSTS_H*/
