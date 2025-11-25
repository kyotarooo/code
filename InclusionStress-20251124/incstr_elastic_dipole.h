#ifndef __INCSTR_ELASTIC_DIPOLE_H_
#define __INCSTR_ELASTIC_DIPOLE_H_

#include "incstr_struct.h"

extern void INCSTR_ElasticDipoleStress(double x[3], int id,
                                       INCLUSION_t *inclusion,
                                       MATERIAL_t *material,
                                       double stress[3][3]);

#endif