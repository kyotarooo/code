#ifndef __DD_DISLOCATION_CORE_H_
#define __DD_DISLOCATION_CORE_H_

#include "dd_struct.h"

extern void DD_DislocationCoreForce(double x[3], double g[3], int coreID,
                                    double *core_x, DD_t *dd, double f[3]);

#endif
