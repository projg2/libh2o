/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_BOUNDARIES_H
#define _H2O_BOUNDARIES_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

enum h2o_region
{
	H2O_REGION_OUT_OF_RANGE,

	H2O_REGION1,
	H2O_REGION2,
	H2O_REGION3,
	H2O_REGION4,
	H2O_REGION5,

	H2O_REGION_MAX
};

enum h2o_region h2o_region_pT(double p, double T);
enum h2o_region h2o_region_ph(double p, double h);
enum h2o_region h2o_region_ps(double p, double s);
enum h2o_region h2o_region_Tx(double T, double x);
enum h2o_region h2o_region_px(double p, double x);
enum h2o_region h2o_region_rhoT(double rho, double T);

double h2o_b23_p_T(double T);
double h2o_b23_T_p(double p);

double h2o_b14_h_s(double s);

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_BOUNDARIES_H*/
