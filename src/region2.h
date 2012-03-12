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

double h2o_region2_v_pT(double p, double T); /* [MPa, K] -> [m³/kg] */
double h2o_region2_u_pT(double p, double T); /* [MPa, K] -> [kJ/kg] */
double h2o_region2_s_pT(double p, double T); /* [MPa, K] -> [kJ/kgK] */
double h2o_region2_h_pT(double p, double T); /* [MPa, K] -> [kJ/kg] */

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_REGION2_H*/
