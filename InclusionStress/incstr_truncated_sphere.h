#ifndef __INCSTR_TRUNCATED_SPHERE_H_
#define __INCSTR_TRUNCATED_SPHERE_H_

#include "incstr_struct.h"

extern void INCSTR_TruncatedSphericalInclusionStress(double x[3], int id,
                                                     INCLUSION_t *inclusion,
                                                     MATERIAL_t *material,
                                                     int *inside,
                                                     double stress[3][3]);

#endif