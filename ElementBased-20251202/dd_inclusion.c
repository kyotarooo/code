#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dd_defs.h"
#include "dd_struct.h"

static double ScalarProduct(double a[3], double b[3]) {
  double c = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    c += a[iDim] * b[iDim];
  }
  return c;
}

static void UnitVector(double a[3]) {
  double r = sqrt(ScalarProduct(a, a));

  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] /= r;
  }
}

static void VectorProduct(double a[3], double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void DD_KeepDislocationsOutsideInclusions(DD_t* dd) {
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x0 = node->x0;
    double* x = node->x;

    if (node->type != MOBILE_NODE) continue;

    for (int i = 0; i < dd->nInclusions; i++) {
      INCLUSION_t* inclusion = &dd->inclusion[i];

      if (inclusion->type == SPHERICAL_INCLUSION) {
        double* xi = inclusion->x;
        double r = inclusion->v[0];
        double r2 = r * r;
        double dx[3];
        double d2;

        for (int iDim = 0; iDim < 3; iDim++) {
          dx[iDim] = x[iDim] - xi[iDim];
        }
        d2 = ScalarProduct(dx, dx);

        // if the dislocation position inside the inclusion
        if (d2 < r2) {
          double t[3];
          double dx0[3];
          double tdx;
          double s;

          for (int iDim = 0; iDim < 3; iDim++) {
            dx[iDim] = x0[iDim] - xi[iDim];
            t[iDim] = x[iDim] - x0[iDim];
            dx0[iDim] = x0[iDim] - xi[iDim];
          }
          d2 = ScalarProduct(dx, dx);
          UnitVector(t);
          tdx = ScalarProduct(t, dx0);
          s = -1.0 * tdx - sqrt(tdx * tdx - (d2 - r2));

          for(int iDim = 0; iDim < 3; iDim++) {
            x[iDim] = x0[iDim] + s * t[iDim];
          }

          break;
        }
      }
    }
  }
}

void DD_SetUpCylindericalInclusionCoordinateSystem(INCLUSION_t *inclusion) {
  double z[3];
  double r[3] = {0.142142124, 0.45644123141, 0.256543634}; // Random vector
  double x[3], y[3];

  inclusion->vector = (double*)malloc(3*sizeof(double));
  inclusion->tensor = (double**)malloc(3*sizeof(double*));

  for(int iDim = 0; iDim < 3; iDim++) {
    z[iDim] = inclusion->v[2 + iDim];
  }
  UnitVector(z);

  // Major axis
  for(int iDim = 0; iDim < 3; iDim++) {
    inclusion->vector[iDim] = z[iDim];
  }

  // Rotation tensor
  VectorProduct(z, r, x);
  UnitVector(x);
  VectorProduct(z, x, y);
  for(int iDim = 0; iDim < 3; iDim++) {
    inclusion->tensor[0][iDim] = x[iDim];
    inclusion->tensor[1][iDim] = y[iDim];
    inclusion->tensor[2][iDim] = z[iDim];
  }
}