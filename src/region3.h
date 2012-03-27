/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_REGION3_H
#define _H2O_REGION3_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

double h2o_region3_p_rhoT(double rho, double T); /* [kg/m³, K] -> [MPa] */
double h2o_region3_u_rhoT(double rho, double T); /* [kg/m³, K] -> [kJ/kg] */
double h2o_region3_s_rhoT(double rho, double T); /* [kg/m³, K] -> [kJ/kgK] */
double h2o_region3_h_rhoT(double rho, double T); /* [kg/m³, K] -> [kJ/kg] */

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_REGION3_H*/
