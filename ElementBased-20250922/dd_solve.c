#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_defs.h"
#include "dd_struct.h"

static void RHSVector(DD_t *dd, double *a) {
  for (int i = 0; i < 3 * dd->nNodes; i++) {
    a[i] = 0.0;
  }

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];
    double *b = element->b;
    for (int iNode = 0; iNode < 2; iNode++) {
      int inodeID = element->nodeID[iNode];
      for (int iDim = 0; iDim < 3; iDim++) {
        a[3 * inodeID + iDim] += b[3 * iNode + iDim];
      }
    }
  }
}

static void DiagonalVector(DD_t *dd, double *a) {
  for (int i = 0; i < 3 * dd->nNodes; i++) {
    a[i] = 0.0;
  }

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];
    double **K = element->K;
    for (int iNode = 0; iNode < 2; iNode++) {
      int inodeID = element->nodeID[iNode];
      for (int iDim = 0; iDim < 3; iDim++) {
        a[3 * inodeID + iDim] += K[3 * iNode + iDim][3 * iNode + iDim];
      }
    }
  }
  for (int i = 0; i < 3 * dd->nNodes; i++) {
    a[i] = 1.0 / sqrt(a[i]);
  }
}

static void BCVector(DD_t *dd, double *a) {
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];
    int type = node->type;

    if (type == IMMOBILE_NODE) {
      a[3 * iNode + 0] = 0.0;
      a[3 * iNode + 1] = 0.0;
      a[3 * iNode + 2] = 0.0;
    } else {
      a[3 * iNode + 0] = 1.0;
      a[3 * iNode + 1] = 1.0;
      a[3 * iNode + 2] = 1.0;
    }
  }
}

static void BoundaryCondition(double *bc, DD_t *dd, double *a) {
  for (int i = 0; i < 3 * dd->nNodes; i++) {
    a[i] *= bc[i];
  }
}

static double ScalarProduct(int n, double *a, double *b) {
  double c = 0.0;

  for (int i = 0; i < n; i++) {
    c += a[i] * b[i];
  }
  return c;
}

static void MapConstraintNodeVector(DD_t *dd, double *a) {
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];
    int type = node->type;

    if (type == CONSTRAINT_NODE || type == FREE_SURFACE_NODE) {
      double *d = node->direction;
      double r = ScalarProduct(3, &a[3 * iNode], d);
      for (int iDim = 0; iDim < 3; iDim++) {
        a[3 * iNode + iDim] = r * d[iDim];
      }
    } else if (type == IMMOBILE_NODE) {
      for (int iDim = 0; iDim < 3; iDim++) {
        a[3 * iNode + iDim] = 0.0;
      }
    }
  }

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];
    double *s = element->slip;

    for (int iNode = 0; iNode < 2; iNode++) {
      int nodeID = element->nodeID[iNode];
      NODE_t *node = &dd->node[nodeID];

      // Remove out-of-slip plane component
      if (node->type == MOBILE_NODE) {
        double r = ScalarProduct(3, s, &a[3 * nodeID]);

        for (int iDim = 0; iDim < 3; iDim++) {
          a[3 * nodeID + iDim] -= r * s[iDim];
        }
      }
    }
  }
}

static void MatrixVectorProduct(double *a, DD_t *dd, double *b) {
  for (int i = 0; i < 3 * dd->nNodes; i++) {
    b[i] = 0.0;
  }

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];
    double **K = element->K;

    for (int iNode = 0; iNode < 2; iNode++) {
      int inodeID = element->nodeID[iNode];
      for (int jNode = 0; jNode < 2; jNode++) {
        int jnodeID = element->nodeID[jNode];
        for (int iDim = 0; iDim < 3; iDim++) {
          for (int jDim = 0; jDim < 3; jDim++) {
            b[3 * inodeID + iDim] +=
                K[3 * iNode + iDim][3 * jNode + jDim] * a[3 * jnodeID + jDim];
          }
        }
      }
    }
  }
}

void DD_SolveEquation(DD_t *dd) {
  int n = 3 * dd->nNodes;
  double *r = (double *)malloc(n * sizeof(double));
  double *d = (double *)malloc(n * sizeof(double));
  double *dr = (double *)malloc(n * sizeof(double));
  double *p = (double *)malloc(n * sizeof(double));
  double *ap = (double *)malloc(n * sizeof(double));
  double *x = (double *)malloc(n * sizeof(double));
  double *b = (double *)malloc(n * sizeof(double));
  double *bc = (double *)malloc(n * sizeof(double));
  double norm = 0.0;
  double lower = 0.0, upper = 0.0;
  double eps = 1.0e-08;

  // Equation of motion is solved using the Conjugate Gradient method with the
  // diagonal scaling preconditioning.
  BCVector(dd, bc);
  DiagonalVector(dd, d);
  RHSVector(dd, b);
  BoundaryCondition(bc, dd, b);

  for (int i = 0; i < n; i++) {
    x[i] = 0.0;
    r[i] = b[i];
    dr[i] = d[i] * r[i];
    p[i] = dr[i];
  }

  lower = ScalarProduct(n, r, dr);
  norm = ScalarProduct(n, r, r);
  if (norm == 0.0) return;

  while (1) {
    double residual = 0.0;
    double c;

    MatrixVectorProduct(p, dd, ap);
    BoundaryCondition(bc, dd, ap);

    upper = lower;
    lower = ScalarProduct(n, p, ap);
    c = upper / lower;

    for (int i = 0; i < n; i++) {
      x[i] += c * p[i];
      r[i] -= c * ap[i];
      dr[i] = d[i] * r[i];
    }

    residual = ScalarProduct(n, r, r);
    if (residual < eps * eps * norm) {
      break;
    }

    lower = ScalarProduct(n, r, dr);
    c = lower / upper;

    for (int i = 0; i < n; i++) {
      p[i] = dr[i] + c * p[i];
    }
  }

  // If the motion of a node is restricted in a specific direction, the solution
  // velocity vector is mapped into the direction.
  MapConstraintNodeVector(dd, x);

  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];
    for (int iDim = 0; iDim < 3; iDim++) {
      node->v[iDim] = x[3 * iNode + iDim];
    }
  }

  free(r);
  free(d);
  free(dr);
  free(p);
  free(ap);
  free(x);
  free(b);
  free(bc);
}
