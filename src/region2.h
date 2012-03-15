/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_REGION2_H
#define _H2O_REGION2_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

enum h2o_region2_subregion
{
	H2O_REGION2_OUT_OF_RANGE,

	H2O_REGION2A,
	H2O_REGION2B,
	H2O_REGION2C,

	H2O_REGION2_MAX
};

double h2o_region2_v_pT(double p, double T); /* [MPa, K] -> [m³/kg] */
double h2o_region2_u_pT(double p, double T); /* [MPa, K] -> [kJ/kg] */
double h2o_region2_s_pT(double p, double T); /* [MPa, K] -> [kJ/kgK] */
double h2o_region2_h_pT(double p, double T); /* [MPa, K] -> [kJ/kg] */

double h2o_region2_b2bc_p_h(double h); /* [kJ/kg] -> [MPa] */
double h2o_region2_b2bc_h_p(double p); /* [MPa] -> [kJ/kg] */

enum h2o_region2_subregion h2o_region2_subregion_ph(double p, double h);

double h2o_region2a_T_ph(double p, double h); /* [MPa, kJ/kg] -> [K] */
double h2o_region2b_T_ph(double p, double h); /* [MPa, kJ/kg] -> [K] */
double h2o_region2c_T_ph(double p, double h); /* [MPa, kJ/kg] -> [K] */

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_REGION2_H*/
