#ifndef __DD_PERIODIC_BC_H_
#define __DD_PERIODIC_BC_H_

#include "dd_struct.h"

extern void DD_PeriodicReturn(double x0[3], double x1[3], BC_t *bc);
extern void DD_PeriodicBoundaryCondition(DD_t *dd);
extern void DD_PeriodicBoundaryConditionPosition(double x[3], DD_t *dd);

#endif
