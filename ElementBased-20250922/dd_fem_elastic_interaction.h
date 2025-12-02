#ifndef __DD_FEM_ELASTIC_INTERACTION_H_
#define __DD_FEM_ELASTIC_INTERACTION_H_

#include "dd_struct.h"
#include "fem_struct.h"

extern void DD_FEM_ElasticInteractionForce(double x[3], double b[3],
                                           double t[3], DD_t *dd, FEM_t *fem,
                                           double f[3]);

#endif
