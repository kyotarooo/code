#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fem_shape.h"
#include "fem_struct.h"

void FEM_NewtonRaphsonMethod(double x[3], double gx[3], FEMELEMENT_t *element,
                             FEM_t *fem) {
  double element_x[20][3];
  double tolerance = 1.0e-12 * 1.0e-12;

  for (int iDim = 0; iDim < 3; iDim++) {
    gx[iDim] = 0.0;
  }

  for (int iNode = 0; iNode < 20; iNode++) {
    int nodeID = element->nodeID[iNode];
    FEMNODE_t *node = &fem->node[nodeID];

    for (int iDim = 0; iDim < 3; iDim++) {
      element_x[iNode][iDim] = node->x[iDim];
    }
  }

  while (1) {
    double shape[20], dshape[20][3];
    double ix[3], dx[3];
    double dx_dg[3][3], dg_dx[3][3];
    double residual = 0.0;
    double det;

    FEM_ShapeFunction3D(gx, shape);

    for (int iDim = 0; iDim < 3; iDim++) {
      double v = 0.0;

      for (int iNode = 0; iNode < 20; iNode++) {
        v += shape[iNode] * element_x[iNode][iDim];
      }
      ix[iDim] = v;
    }

    for (int iDim = 0; iDim < 3; iDim++) {
      dx[iDim] = x[iDim] - ix[iDim];
      residual += dx[iDim] * dx[iDim];
    }

    if (residual < tolerance) return;

    // update
    FEM_DShapeFunction3D(gx, dshape);
    for (int iDim = 0; iDim < 3; iDim++) {
      for (int jDim = 0; jDim < 3; jDim++) {
        double v = 0.0;

        for (int iNode = 0; iNode < 20; iNode++) {
          v += dshape[iNode][jDim] * element_x[iNode][iDim];
        }
        dx_dg[iDim][jDim] = v;
      }
    }
    det =
        1.0 /
        (dx_dg[0][0] * (dx_dg[1][1] * dx_dg[2][2] - dx_dg[1][2] * dx_dg[2][1]) +
         dx_dg[0][1] * (dx_dg[1][2] * dx_dg[2][0] - dx_dg[1][0] * dx_dg[2][2]) +
         dx_dg[0][2] * (dx_dg[1][0] * dx_dg[2][1] - dx_dg[1][1] * dx_dg[2][0]));
    dg_dx[0][0] = (dx_dg[1][1] * dx_dg[2][2] - dx_dg[1][2] * dx_dg[2][1]) * det;
    dg_dx[0][1] = (dx_dg[0][2] * dx_dg[2][1] - dx_dg[0][1] * dx_dg[2][2]) * det;
    dg_dx[0][2] = (dx_dg[0][1] * dx_dg[1][2] - dx_dg[0][2] * dx_dg[1][1]) * det;
    dg_dx[1][0] = (dx_dg[1][2] * dx_dg[2][0] - dx_dg[1][0] * dx_dg[2][2]) * det;
    dg_dx[1][1] = (dx_dg[0][0] * dx_dg[2][2] - dx_dg[0][2] * dx_dg[2][0]) * det;
    dg_dx[1][2] = (dx_dg[0][2] * dx_dg[1][0] - dx_dg[0][0] * dx_dg[1][2]) * det;
    dg_dx[2][0] = (dx_dg[1][0] * dx_dg[2][1] - dx_dg[1][1] * dx_dg[2][0]) * det;
    dg_dx[2][1] = (dx_dg[0][1] * dx_dg[2][0] - dx_dg[0][0] * dx_dg[2][1]) * det;
    dg_dx[2][2] = (dx_dg[0][0] * dx_dg[1][1] - dx_dg[0][1] * dx_dg[1][0]) * det;

    for (int iDim = 0; iDim < 3; iDim++) {
      double v = 0.0;

      for (int jDim = 0; jDim < 3; jDim++) {
        v += dg_dx[iDim][jDim] * dx[jDim];
      }
      gx[iDim] += v;
    }
  }
}
