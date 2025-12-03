#include <math.h>
#include <stdio.h>

#include "dd_struct.h"

static int InsideInclusion(double x[3], double g[3], int coreID, DD_t *dd) {
  int nDislocations = dd->core.nDislocations;
  double *templete_x = dd->core.property[coreID].x;
  double width = templete_x[nDislocations - 1] - templete_x[0];

  for (int i = 0; i < dd->nInclusions; i++) {
    INCLUSION_t *inclusion = &dd->inclusion[i];
  }
  return 0;
}

void DD_DislocationCoreForce(double x[3], double g[3], int coreID,
                             double *core_x, DD_t *dd, double f[3]) {
  int nDislocations = dd->core.nDislocations;
  double be[1000];
  double bs[1000];
  double fmatrix[1000], finclusion[1000];

  if (nDislocations == 0) {
    for (int iDim = 0; iDim < 3; iDim++) {
      f[iDim] = 0.0;
    }
    return;
  }

  if (!InsideInclusion(x, g, coreID, dd)) {
    double residual = 0.0;
    for (int i = 0; i < nDislocations; i++) {
      residual += dd->core.property[coreID].fmatrix[i];
    }
    for (int iDim = 0; iDim < 3; iDim++) {
      f[iDim] = residual * g[iDim];
    }
    return;
  }

  for (int i = 0; i < nDislocations; i++) {
    be[i] = dd->core.property[coreID].bedge[i];
    bs[i] = dd->core.property[coreID].bscrew[i];
    fmatrix[i] = dd->core.property[coreID].fmatrix[i];
    finclusion[i] = dd->core.property[coreID].finclusion[i];
  }
}
