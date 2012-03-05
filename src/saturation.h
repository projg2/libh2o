/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_SATURATION_H
#define _H2O_SATURATION_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

double h2o_saturation_p_T(double T); /* p [MPa] = f(T [K]) */
double h2o_saturation_T_p(double p); /* T [K] = f(p [MPa]) */

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_SATURATION_H*/
