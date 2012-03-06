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

static inline double pow2(double arg)
{
	return arg * arg;
}

static inline double pow4(double arg)
{
	return pow2(pow2(arg));
}

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_XMATH_H*/
