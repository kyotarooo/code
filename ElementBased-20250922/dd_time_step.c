#include <math.h>
#include <stdio.h>

#include "dd_struct.h"

double DD_AdjustTimeIncrement(DD_t *dd) {
  double vmax = 0.0;
  double dtmin = dd->step.dtmin;
  double dtmax = dd->step.dtmax;
  double dt;

  // Find the maximum velocity
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];
    double v = 0.0;
    for (int iDim = 0; iDim < 3; iDim++) {
      v += node->v[iDim] * node->v[iDim];
    }
    vmax = fmax(vmax, v);
  }
  vmax = sqrt(vmax);

  // Control the maximum displacement per time step
  dt = dd->step.umax / vmax;

  // The time increment must be between the minimum and the maximum
  if (dtmin > 0.0 && dt < dtmin)
    dt = dtmin;
  else if (dtmax > 0.0 && dt > dtmax)
    dt = dtmax;

  return dt;
}