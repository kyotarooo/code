#include <stdio.h>
#include <stdlib.h>

#include "dd_struct.h"

void DD_PeriodicReturn(double x0[3], double x1[3], BC_t *bc) {
  for (int iDim = 0; iDim < 3; iDim++) {
    double size = bc->size[iDim];

    if (bc->periodic[iDim]) {
      double dx = x1[iDim] - x0[iDim];

      if (dx < -0.5 * size) {
        x1[iDim] += size;
      } else if (dx > 0.5 * size) {
        x1[iDim] -= size;
      }
    }
  }
}

void DD_PeriodicBoundaryCondition(DD_t *dd) {
  for (int iDim = 0; iDim < 3; iDim++) {
    if (dd->bc.periodic[iDim]) {
      double size = dd->bc.size[iDim];

      for (int iNode = 0; iNode < dd->nNodes; iNode++) {
        NODE_t *node = &dd->node[iNode];
        double x = node->x[iDim];

        if (x < 0.0)
          node->x[iDim] += size;
        else if (x > size)
          node->x[iDim] -= size;
      }
    }
  }
}

void DD_PeriodicBoundaryConditionPosition(double x[3], DD_t *dd) {
  for (int iDim = 0; iDim < 3; iDim++) {
    if (dd->bc.periodic[iDim]) {
      double xi = x[iDim];
      double size = dd->bc.size[iDim];

      if (xi < 0.0)
        x[iDim] += size;
      else if (xi > size)
        x[iDim] -= size;
    }
  }
}
