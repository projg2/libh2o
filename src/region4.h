/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_REGION4_H
#define _H2O_REGION4_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

double h2o_region4_v_Tx(double T, double x); /* [K, 0..1] -> [m³/kg] */
double h2o_region4_u_Tx(double T, double x); /* [K, 0..1] -> [kJ/kg] */
double h2o_region4_s_Tx(double T, double x); /* [K, 0..1] -> [kJ/kgK] */
double h2o_region4_h_Tx(double T, double x); /* [K, 0..1] -> [kJ/kg] */

double h2o_region4_x_Ts(double T, double s); /* [K, kJ/kgK] -> [0..1] */
double h2o_region4_x_Th(double T, double h); /* [K, kJ/kg] -> [0..1] */

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_REGION4_H*/
