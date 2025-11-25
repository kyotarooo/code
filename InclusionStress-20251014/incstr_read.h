#ifndef __INCSTR_READ_H_
#define __INCSTR_READ_H_

#include "incstr_struct.h"

extern void INCSTR_ReadInclusions(INCLUSION_t *inclusion, char *directory);

extern void INCSTR_ReadGrids(GRID_CELL_t *grid_cell, MATERIAL_t *material,
                             char *directory);

extern void INCSTR_ReadMaterial(MATERIAL_t *material, char *directory);

#endif