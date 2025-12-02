#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_constraint.h"
#include "dd_defs.h"
#include "dd_fem_velocity.h"
#include "dd_node_handle.h"
#include "dd_struct.h"
#include "dd_time_step.h"
#include "fem_struct.h"

static void PeriodicBoundaryCondition(DD_t *dd) {
  for (int iDim = 0; iDim < 3; iDim++) {
    if (dd->bc.periodic[iDim]) {
      double size = dd->bc.size[iDim];

      for (int iNode = 0; iNode < dd->nNodes; iNode++) {
        NODE_t *node = &dd->node[iNode];
        double x = node->x[iDim];

        if (x < 0.0) {
          node->x[iDim] += size;
          node->v2[iDim] += size;
        } else if (x > size) {
          node->x[iDim] -= size;
          node->v2[iDim] -= size;
        }
      }
    }
  }
}

static void UpdateInitials(DD_t *dd) {
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];
    for (int iDim = 0; iDim < 3; iDim++) {
      node->v0[iDim] = node->x[iDim];
      node->v1[iDim] = node->v[iDim];
    }
  }
}

static void TrapezoidalIntegral(double dt, DD_t *dd, FEM_t *fem) {
  double tolerance = dd->bc.element_length[1] * TRAPEZOIDAL_EPS;

  // v0: Position at the current configuration
  // v1: Velocity at the current configuration
  // v2: Tentative position for trapezoidal integral
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];
    for (int iDim = 0; iDim < 3; iDim++) {
      node->x[iDim] = node->v0[iDim] + dt * node->v1[iDim];
      node->v2[iDim] = node->v0[iDim] + 0.5 * dt * node->v1[iDim];
    }
  }
  DD_ConstrainNodes(dd);
  PeriodicBoundaryCondition(dd);

  for (int i = 0; i < NMAX_TRAPEZOIDAL_ITERATIONS; i++) {
    double e = 0.0;

    DD_FEM_DislocationVelocity(dd, fem);
    for (int iNode = 0; iNode < dd->nNodes; iNode++) {
      NODE_t *node = &dd->node[iNode];
      double ei = 0.0;
      for (int iDim = 0; iDim < 3; iDim++) {
        double x = node->v2[iDim] + 0.5 * dt * node->v[iDim];
        double dx = node->x[iDim] - x;

        node->x[iDim] = x;
        ei += dx * dx;
      }
      if (node->type == MOBILE_NODE) {
        e = fmax(e, ei);  // Error is controlled only for mobile nodes
      }
    }
    DD_ConstrainNodes(dd);
    PeriodicBoundaryCondition(dd);

    if (e < tolerance * tolerance) {
      return;  // The solution could be obtained for this time interval
    }
  }

  // Subdivide time step size
  dt /= (double)N_TRAPEZOIDAL_SUBDIVISIONS;
  for (int i = 0; i < N_TRAPEZOIDAL_SUBDIVISIONS; i++) {
    TrapezoidalIntegral(dt, dd, fem);
    UpdateInitials(dd);
  }
}

void DD_FEM_TrapezoidalIntegral(DD_t *dd, FEM_t *fem) {
  double dt;

  DD_BackupPosition(dd);
  DD_FindElementsAroundNode(dd);

  DD_FEM_DislocationVelocity(dd, fem);
  dt = DD_AdjustTimeIncrement(dd);

  UpdateInitials(dd);
  TrapezoidalIntegral(dt, dd, fem);
  dd->step.dt = dt;
  dd->step.t += dt;
}