#include <math.h>
#include <stdio.h>

#include "dd_defs.h"
#include "dd_struct.h"

void DD_ConstrainNodes(DD_t *dd) {
  double *size = dd->bc.size;

  return;
  
  // To keep the nodal positions within the simulation volume
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];
    double *x0 = node->x0;
    double *x = node->x;
    double dx[3];

    if (node->type != MOBILE_NODE) continue;

    for (int iDim = 0; iDim < 3; iDim++) {
      dx[iDim] = x[iDim] - x0[iDim];
    }

    for (int iDim = 0; iDim < 3; iDim++) {
      if (x[iDim] < 0.0) {
        double t = (0.0 - x0[iDim]) / dx[iDim];

        for(int iDim = 0; iDim < 3; iDim++) {
          x[iDim] = x0[iDim] + t * dx[iDim];
        }
      } else if (x[iDim] > size[iDim]) {
        double t = (size[iDim] - x0[iDim]) / dx[iDim];

        for(int iDim = 0; iDim < 3; iDim++) {
          x[iDim] = x0[iDim] + t * dx[iDim];
        }
      }
    }
  }
}