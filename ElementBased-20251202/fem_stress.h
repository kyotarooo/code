#ifndef __FEM_STRESS_H_
#define __FEM_STRESS_H_

#include "fem_struct.h"

extern void FEM_ComputeStress(FEM_t *fem);
extern void FEM_CorrectionStress(double x[3], FEM_t *fem, double stress[3][3]);

#endif
