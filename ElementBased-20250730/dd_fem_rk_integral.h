#ifndef __DD_FEM_RK_INTEGRAL_H_
#define __DD_FEM_RK_INTEGRAL_H_

#include "dd_struct.h"
#include "fem_struct.h"

extern void DD_FEM_RungeKutta1stIntegral(DD_t *dd, FEM_t *fem);
extern void DD_FEM_RungeKutta2ndIntegral(DD_t *dd, FEM_t *fem);
extern void DD_FEM_ungeKutta3rdIntegral(DD_t *dd, FEM_t *fem);
extern void DD_FEM_RungeKutta4thIntegral(DD_t *dd, FEM_t *fem);

#endif
