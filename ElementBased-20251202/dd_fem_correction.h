#ifndef __DD_FEM_CORRECTION_H_
#define __DD_FEM_CORRECTION_H_

#include "dd_struct.h"
#include "fem_struct.h"

extern void DD_FEM_SurfaceTractionGridStress(DD_t *dd, FEM_t *fem);
extern void DD_FEM_CorrectionField(DD_t *dd, FEM_t *fem);

#endif
