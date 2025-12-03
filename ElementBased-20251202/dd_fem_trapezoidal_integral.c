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

static void PeriodicBoundaryCondition(DD_t* dd) {
  for (int iDim = 0; iDim < 3; iDim++) {
    if (dd->bc.periodic[iDim]) {
      double size = dd->bc.size[iDim];

      for (int iNode = 0; iNode < dd->nNodes; iNode++) {
        NODE_t* node = &dd->node[iNode];
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

static void UpdateInitials(DD_t* dd) {
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* v = node->v;
    double* v0 = node->v0;
    double* v1 = node->v1;
    for (int iDim = 0; iDim < 3; iDim++) {
      v0[iDim] = x[iDim];
      v1[iDim] = v[iDim];
    }
  }
}

static void TrapezoidalIntegral(double dt, DD_t* dd, FEM_t* fem) {
  double tolerance = dd->bc.element_length[1] * TRAPEZOIDAL_EPS;
  // This parameter controls the time integral scheme:
  // 0.5: Semi-implicit (Trapezoidal rule), 1.0: Implicit (Backward Euler)
  double theta = dd->step.theta;

  // v0: Position at the current configuration
  // v1: Velocity at the current configuration
  // v2: Tentative position for trapezoidal integral
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];
    double* x = node->x;
    double* v0 = node->v0;
    double* v1 = node->v1;
    double* v2 = node->v2;
    for (int iDim = 0; iDim < 3; iDim++) {
      x[iDim] = v0[iDim] + dt * v1[iDim];
      v2[iDim] = v0[iDim] + (1.0 - theta) * dt * v1[iDim];
    }
  }
  DD_ConstrainNodes(dd);
  PeriodicBoundaryCondition(dd);

  for (int i = 0; i < NMAX_TRAPEZOIDAL_ITERATIONS; i++) {
    double e = 0.0;

    DD_FEM_DislocationVelocity(dd, fem);
    for (int iNode = 0; iNode < dd->nNodes; iNode++) {
      NODE_t* node = &dd->node[iNode];
      double* x = node->x;
      double* v = node->v;
      double* v1 = node->v1;
      double* v2 = node->v2;
      double* v3 = node->v3;
      for (int iDim = 0; iDim < 3; iDim++) {
        v3[iDim] = x[iDim];
        x[iDim] = v2[iDim] + theta * dt * v[iDim];
        v[iDim] = (1.0 - theta) * v1[iDim] + theta * v[iDim];
      }
    }
    DD_ConstrainNodes(dd);
    for (int iNode = 0; iNode < dd->nNodes; iNode++) {
      NODE_t* node = &dd->node[iNode];
      double* x = node->x;
      double* v3 = node->v3;
      double ei = 0.0;
      for (int iDim = 0; iDim < 3; iDim++) {
        double dx = v3[iDim] - x[iDim];

        ei += dx * dx;
        if (node->type == MOBILE_NODE) {
          e = fmax(e, ei);  // Error is controlled only for mobile nodes
        }
      }
    }
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

void DD_FEM_TrapezoidalIntegral(DD_t* dd, FEM_t* fem) {
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