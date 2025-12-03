#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_defs.h"
#include "dd_dislocation_core.h"
#include "dd_element_handle.h"
#include "dd_fem_elastic_interaction.h"
#include "dd_periodic_bc.h"
#include "dd_shape.h"
#include "dd_solve.h"
#include "dd_struct.h"

static void VectorProduct(double a[3], double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static void InitializeTensor(int ni, int nj, double **a) {
  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      a[i][j] = 0.0;
    }
  }
}

static void InitializeVector(int n, double *a) {
  for (int i = 0; i < n; i++) {
    a[i] = 0.0;
  }
}

static void SubstituteVector(double *a, double b[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] = a[iDim];
  }
}

static void IntegralPointPosition(double x0[3], double x1[3], double shape[2],
                                  double x[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    x[iDim] = shape[0] * x0[iDim] + shape[1] * x1[iDim];
  }
}

static void AddVector(double a[3], double b[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] += a[iDim];
  }
}

static void ProjectVector(double a[3], double b[3]) {
  double r = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    r += a[iDim] * b[iDim];
  }
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] = r * a[iDim];
  }
}

static void DislocationElementIntegral(int elementID, DD_t *dd, FEM_t *fem) {
  ELEMENT_t *element = &dd->element[elementID];
  double *element_b = element->b;
  double **element_K = element->K;
  int coreID = element->coreID;
  int nodeID0 = element->nodeID[0];
  int nodeID1 = element->nodeID[1];
  NODE_t *node0 = &dd->node[nodeID0];
  NODE_t *node1 = &dd->node[nodeID1];
  double x0[3], x1[3];
  double b[3], s[3], t[3], g[3];
  double jacobian = 0.5 * element->length;  // -1 < u < 1 -> du = 2
  double c = 1.0 / dd->material.mobility;
  double gauss_u[4] = {
      -0.861136311594052572537805190223, -0.339981043584856312822495283399,
      0.339981043584856312822495283399, 0.861136311594052572537805190223};
  double gauss_w[4] = {
      0.347854845137453849712727560473, 0.652145154862546205798423670785,
      0.652145154862546205798423670785, 0.347854845137453849712727560473};

  // Initialize
  InitializeTensor(6, 6, element->K);
  InitializeVector(6, element->b);

  // Element information
  SubstituteVector(node0->x, x0);
  SubstituteVector(node1->x, x1);
  SubstituteVector(element->burgers, b);
  SubstituteVector(element->slip, s);
  SubstituteVector(element->tangent, t);
  DD_PeriodicReturn(x0, x1, &dd->bc);

  // Glide direction
  VectorProduct(s, t, g);

  for (int iIntegral = 0; iIntegral < N_DD_INTEGRALS; iIntegral++) {
    double u = gauss_u[iIntegral];
    double w = gauss_w[iIntegral];
    double shape[2], x[3];
    double felastic[3] = {0.0, 0.0, 0.0};
    double fcore[3] = {0.0, 0.0, 0.0};
    double f[3] = {0.0, 0.0, 0.0};
    double jw, cjw;
    double *core_x = element->integral[iIntegral].core_x;

    // Integral point position
    DD_ShapeFunction1D(u, shape);
    IntegralPointPosition(x0, x1, shape, x);

    // Elastic interaction force: felastic
    // Applied, dislocation, grid and correction stress
    DD_FEM_ElasticInteractionForce(x, b, t, dd, fem, felastic);

    // Dislocation core force: fcore
    // Virtual dislocation core model
    DD_DislocationCoreForce(x, g, coreID, core_x, dd, fcore);

    // Total force and project force vector to the glide direction
    AddVector(felastic, f);
    AddVector(fcore, f);
    ProjectVector(g, f);

    // RHS vector
    jw = jacobian * w;
    for (int iNode = 0; iNode < 2; iNode++) {
      double sjw = shape[iNode] * jw;
      for (int iDim = 0; iDim < 3; iDim++) {
        element_b[3 * iNode + iDim] += f[iDim] * sjw;
      }
    }

    // Stiffness (resistivity) matrix
    cjw = c * jw;
    for (int iNode = 0; iNode < 2; iNode++) {
      double csjw = shape[iNode] * cjw;
      for (int jNode = 0; jNode < 2; jNode++) {
        double cssjw = shape[jNode] * csjw;
        for (int iDim = 0; iDim < 3; iDim++) {
          element_K[3 * iNode + iDim][3 * jNode + iDim] += cssjw;
        }
      }
    }
  }
}

void DD_FEM_DislocationVelocity(DD_t *dd, FEM_t *fem) {
  // Update all tangent vectors
  DD_TangentVectors(dd);

#pragma omp parallel for
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    DislocationElementIntegral(iElement, dd, fem);
  }
  DD_SolveEquation(dd);
}
