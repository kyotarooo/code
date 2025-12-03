#ifndef __INCSTR_DILATATION_CENTER_H_
#define __INCSTR_DILATATION_CENTER_H_

#include "incstr_struct.h"

extern void INCSTR_DilatationCenterStress(double x[3], int id,
                                       INCLUSION_t *inclusion,
                                       MATERIAL_t *material,
                                       double stress[3][3]);

#endif