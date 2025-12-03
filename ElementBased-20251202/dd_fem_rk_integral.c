#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_constraint.h"
#include "dd_fem_velocity.h"
#include "dd_node_handle.h"
#include "dd_periodic_bc.h"
#include "dd_struct.h"
#include "dd_time_step.h"
#include "fem_struct.h"

void DD_FEM_RungeKutta1stIntegral(DD_t* dd, FEM_t* fem) {
  double dt;

  DD_BackupPosition(dd);
  DD_FindElementsAroundNode(dd);

  // 1st step
  DD_FEM_DislocationVelocity(dd, fem);
  dd->step.dt = dt = DD_AdjustTimeIncrement(dd);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v = node->v;
    for (int iDim = 0; iDim < 3; iDim++) {
      x[iDim] = x0[iDim] + v[iDim] * dt;
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);

  dd->step.t += dt;
}

void DD_FEM_RungeKutta2ndIntegral(DD_t* dd, FEM_t* fem) {
  double dt;
  double c0 = 0.5;
  double c1 = 0.5;

  DD_BackupPosition(dd);
  DD_FindElementsAroundNode(dd);

  // 1st step
  DD_FEM_DislocationVelocity(dd, fem);
  dd->step.dt = dt = DD_AdjustTimeIncrement(dd);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v = node->v;
    double* v0 = node->v0;
    for (int iDim = 0; iDim < 3; iDim++) {
      x[iDim] = x0[iDim] + v[iDim] * dt;
      v0[iDim] = v[iDim];
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);
  // 2nd step
  DD_FEM_DislocationVelocity(dd, fem);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v0 = node->v0;
    double* v = node->v;
    for (int iDim = 0; iDim < 3; iDim++) {
      double vi = c0 * v0[iDim] + c1 * v[iDim];

      x[iDim] = x0[iDim] + vi * dt;
      v[iDim] = vi;
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);

  dd->step.t += dt;
}

void DD_FEM_RungeKutta3rdIntegral(DD_t* dd, FEM_t* fem) {
  double dt;
  double dth, dt2;
  double c0 = 1.0 / 6.0;
  double c1 = 2.0 / 3.0;
  double c2 = 1.0 / 6.0;

  DD_BackupPosition(dd);
  DD_FindElementsAroundNode(dd);

  // 1st step
  DD_FEM_DislocationVelocity(dd, fem);
  dd->step.dt = dt = DD_AdjustTimeIncrement(dd);
  dth = 0.5 * dt;
  dt2 = 2.0 * dt;
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v = node->v;
    double* v0 = node->v0;
    for (int iDim = 0; iDim < 3; iDim++) {
      x[iDim] = x0[iDim] + v[iDim] * dth;
      v0[iDim] = v[iDim];
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);
  // 2nd step
  DD_FEM_DislocationVelocity(dd, fem);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v = node->v;
    double* v0 = node->v0;
    double* v1 = node->v1;
    for (int iDim = 0; iDim < 3; iDim++) {
      x[iDim] = x0[iDim] - v0[iDim] * dt + v[iDim] * dt2;
      v1[iDim] = v[iDim];
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);
  // 3rd step
  DD_FEM_DislocationVelocity(dd, fem);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v = node->v;
    double* v0 = node->v0;
    double* v1 = node->v1;
    for (int iDim = 0; iDim < 3; iDim++) {
      double vi = c0 * v0[iDim] + c1 * v1[iDim] + c2 * v[iDim];

      x[iDim] = x0[iDim] + vi * dt;
      v[iDim] = vi;
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);

  dd->step.t += dt;
}

void DD_FEM_RungeKutta4thIntegral(DD_t* dd, FEM_t* fem) {
  double dt;
  double dth;
  double c0 = 1.0 / 6.0;
  double c1 = 1.0 / 3.0;
  double c2 = 1.0 / 3.0;
  double c3 = 1.0 / 6.0;

  DD_BackupPosition(dd);
  DD_FindElementsAroundNode(dd);

  // 1st step
  DD_FEM_DislocationVelocity(dd, fem);
  dd->step.dt = dt = DD_AdjustTimeIncrement(dd);
  dth = 0.5 * dt;
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v = node->v;
    double* v0 = node->v0;
    for (int iDim = 0; iDim < 3; iDim++) {
      x[iDim] = x0[iDim] + v[iDim] * dth;
      v0[iDim] = v[iDim];
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);
  // 2nd step
  DD_FEM_DislocationVelocity(dd, fem);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v = node->v;
    double* v1 = node->v1;
    for (int iDim = 0; iDim < 3; iDim++) {
      x[iDim] = x0[iDim] + v[iDim] * dth;
      v1[iDim] = v[iDim];
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);
  // 3rd step
  DD_FEM_DislocationVelocity(dd, fem);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v = node->v;
    double* v2 = node->v2;
    for (int iDim = 0; iDim < 3; iDim++) {
      x[iDim] = x0[iDim] + v[iDim] * dt;
      v2[iDim] = v[iDim];
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);
  // 4th step
  DD_FEM_DislocationVelocity(dd, fem);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* x0 = node->x0;
    double* v = node->v;
    double* v0 = node->v0;
    double* v1 = node->v1;
    double* v2 = node->v2;
    for (int iDim = 0; iDim < 3; iDim++) {
      double vi = c0 * v0[iDim] + c1 * v1[iDim] + c2 * v2[iDim] + c3 * v[iDim];

      x[iDim] = x0[iDim] + vi * dt;
      v[iDim] = vi;
    }
  }
  DD_ConstrainNodes(dd);
  DD_PeriodicBoundaryCondition(dd);

  dd->step.t += dt;
}
