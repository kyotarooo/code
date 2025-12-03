#ifndef __FEM_SOLVE_H_
#define __FEM_SOLVE_H_

#include "fem_struct.h"

extern void FEM_Factorization(FEM_t *fem);
extern void FEM_ComputeSolution(FEM_t *fem);

#endif
