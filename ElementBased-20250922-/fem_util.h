#ifndef __FEM_UTIL_H_
#define __FEM_UTIL_H_

#include "fem_struct.h"

extern void FEM_NewtonRaphsonMethod(double x[3], double gx[3],
                                    FEMELEMENT_t *element, FEM_t *fem);

#endif
