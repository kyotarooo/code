#ifndef __FEM_WRITE_H_
#define __FEM_WRITE_H_

#include "fem_struct.h"

extern void FEM_WriteNodes(FEM_t *fem, char *directory);
extern void FEM_WriteMesh(FEM_t *fem, char *directory);
extern void FEM_WriteStress(int stressID, FEM_t *fem, char *directory);

#endif
