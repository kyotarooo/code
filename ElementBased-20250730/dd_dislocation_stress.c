#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_defs.h"
#include "dd_periodic_bc.h"
#include "dd_struct.h"

static void InitializeTensor(double a[3][3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      a[iDim][jDim] = 0.0;
    }
  }
}

static void SubstituteVector(double *a, double b[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] = a[iDim];
  }
}

static double ScalarProduct(double a[3], double b[3]) {
  double c = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    c += a[iDim] * b[iDim];
  }
  return c;
}

static void VectorProduct(double a[3], double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static void TensorProduct(double a[3], double b[3], double c[3][3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = iDim; jDim < 3; jDim++) {
      c[iDim][jDim] = a[iDim] * b[jDim];
    }
  }
}

static void SymmetricTensorProduct(double a[3], double b[3], double c[3][3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = iDim; jDim < 3; jDim++) {
      c[iDim][jDim] = a[iDim] * b[jDim] + a[jDim] * b[iDim];
    }
  }
}

static void SymmetricTensor(double a[3][3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = iDim + 1; jDim < 3; jDim++) {
      a[jDim][iDim] = a[iDim][jDim];
    }
  }
}

static void DistanceVector(double x[3], double xp[3], double t[3],
                           double d[3]) {
  double R[3];
  double Rdt;

  for (int iDim = 0; iDim < 3; iDim++) {
    R[iDim] = x[iDim] - xp[iDim];
  }
  Rdt = ScalarProduct(R, t);
  for (int iDim = 0; iDim < 3; iDim++) {
    d[iDim] = R[iDim] - Rdt * t[iDim];
  }
}

static double SquaredVectorLength(double a[3]) {
  double length2 = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    length2 += a[iDim] * a[iDim];
  }
  return length2;
}

static void CoefficientsForFiniteDislocationElement(double a2d2, double x[3],
                                                    double xp[3], double t[3],
                                                    double c[5]) {
  double a2d2inv = 1.0 / a2d2;
  double R[3];
  double Rdt;
  double Rainv, Ra3inv;
  double s;

  for (int iDim = 0; iDim < 3; iDim++) {
    R[iDim] = x[iDim] - xp[iDim];
  }
  Rdt = ScalarProduct(R, t);
  s = -1.0 * Rdt;

  Rainv = 1.0 / sqrt(a2d2 + Rdt * Rdt);
  Ra3inv = Rainv * Rainv * Rainv;

  c[0] = s * Rainv * a2d2inv;
  c[1] = -1.0 * Rainv;
  c[2] = (2.0 * c[0] + s * Ra3inv) * a2d2inv;
  c[3] = -1.0 * Ra3inv;
  c[4] = c[0] - s * Ra3inv;
}

static void FiniteDislocationElementStress(int elementID, double x[3], DD_t *dd,
                                           double stress[3][3]) {
  ELEMENT_t *element = &dd->element[elementID];
  int nodeID0 = element->nodeID[0];
  int nodeID1 = element->nodeID[1];
  NODE_t *node0 = &dd->node[nodeID0];
  NODE_t *node1 = &dd->node[nodeID1];
  double g = dd->material.g;
  double v = dd->material.v;
  double a = dd->material.a_core;
  double a2 = a * a;
  double g4p = g / (4.0 * M_PI);
  double g4pv = g4p / (1.0 - v);
  double g8p = g4p / 2.0;
  double x0[3], x1[3];
  double t[3], b[3];
  double d[3], d2 = 0.0;
  double c0[5], c1[5], c[5];
  double dxb[3], dxbdt;
  double txb[3];
  double tct[3][3], dcd[3][3];
  double tctxb_sym[3][3];
  double tcd_sym[3][3];
  double tcdxb_sym[3][3];
  double dctxb_sym[3][3];
  double I[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  // Position of edge nodes
  SubstituteVector(node0->x, x0);
  SubstituteVector(node1->x, x1);
  DD_PeriodicReturn(x, x0, &dd->bc);
  DD_PeriodicReturn(x0, x1, &dd->bc);

  // Element information
  SubstituteVector(element->tangent, t);
  SubstituteVector(element->burgers, b);

  // Distance vector from the element
  DistanceVector(x, x0, t, d);
  d2 = SquaredVectorLength(d);

  // Coefficients
  CoefficientsForFiniteDislocationElement(a2 + d2, x, x0, t, c0);
  CoefficientsForFiniteDislocationElement(a2 + d2, x, x1, t, c1);
  for (int i = 0; i < 5; i++) {
    c[i] = c1[i] - c0[i];
  }

  VectorProduct(d, b, dxb);
  VectorProduct(t, b, txb);
  TensorProduct(t, t, tct);
  TensorProduct(d, d, dcd);
  SymmetricTensorProduct(t, txb, tctxb_sym);
  SymmetricTensorProduct(t, d, tcd_sym);
  SymmetricTensorProduct(t, dxb, tcdxb_sym);
  SymmetricTensorProduct(d, txb, dctxb_sym);
  dxbdt = ScalarProduct(dxb, t);

  // Loop for symmetric
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = iDim; jDim < 3; jDim++) {
      stress[iDim][jDim] +=
          c[0] * (g4pv * dxbdt * I[iDim][jDim] + g4pv * dctxb_sym[iDim][jDim] -
                  g4p * tcdxb_sym[iDim][jDim]) +
          c[1] * (-1.0 * g4pv * v * tctxb_sym[iDim][jDim]) +
          c[2] * (g4pv * dxbdt * dcd[iDim][jDim] +
                  g4pv * dxbdt * a2 * I[iDim][jDim] -
                  g8p * a2 * tcdxb_sym[iDim][jDim]) +
          c[3] * (g8p * a2 * tctxb_sym[iDim][jDim] -
                  g4pv * dxbdt * tcd_sym[iDim][jDim]) +
          c[4] * (g4pv * dxbdt * tct[iDim][jDim]);
    }
  }
}

#ifdef _NOTERMINATION_
static void OppositeVector(double a[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] *= -1.0;
  }
}

static void CoefficientsForSemiInfiniteDislocationElement(
    double a2d2, double x[3], double xp[3], double t[3], double c[5]) {
  double a2d2inv = 1.0 / a2d2;
  double R[3];
  double Rdt;
  double Rainv, Ra3inv;
  double s;

  for (int iDim = 0; iDim < 3; iDim++) {
    R[iDim] = x[iDim] - xp[iDim];
  }
  Rdt = ScalarProduct(R, t);
  s = -1.0 * Rdt;

  Rainv = 1.0 / sqrt(a2d2 + Rdt * Rdt);
  Ra3inv = Rainv * Rainv * Rainv;

  c[0] = (1.0 - s * Rainv) * a2d2inv;
  c[1] = Rainv;
  c[2] = (2.0 * c[0] - s * Ra3inv) * a2d2inv;
  c[3] = Ra3inv;
  c[4] = c[0] + s * Ra3inv;
}

static void SemiInfiniteDislocationElementStress(int elementID, double x[3],
                                                 DD_t *dd,
                                                 double stress[3][3]) {
  ELEMENT_t *element = &dd->element[elementID];
  int nodeID0 = element->nodeID[0];
  int nodeID1 = element->nodeID[1];
  NODE_t *node0 = &dd->node[nodeID0];
  NODE_t *node1 = &dd->node[nodeID1];
  double g = dd->material.g;
  double v = dd->material.v;
  double a = dd->material.a_core;
  double g4p = g / (4.0 * M_PI);
  double g4pv = g4p / (1.0 - v);
  double g8p = g4p / 2.0;
  double a2 = a * a;
  double x0[3], x1[3];
  double t[3], b[3];
  double d[3], d2 = 0.0;
  double c[5];
  double dxb[3], dxbdt;
  double txb[3];
  double tct[3][3], dcd[3][3];
  double tctxb_sym[3][3];
  double tcd_sym[3][3];
  double tcdxb_sym[3][3];
  double dctxb_sym[3][3];
  double I[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  // Position of edge nodes
  SubstituteVector(node0->x, x0);
  SubstituteVector(node1->x, x1);
  DD_PeriodicReturn(x, x0, &dd->bc);
  DD_PeriodicReturn(x0, x1, &dd->bc);

  // Element information
  SubstituteVector(element->tangent, t);
  SubstituteVector(element->burgers, b);

  // Distance vector from the element
  DistanceVector(x, x0, t, d);
  d2 = SquaredVectorLength(d);

  // Coefficients
  if (node0->nElements == 1) {
    OppositeVector(t);
    OppositeVector(b);
    CoefficientsForSemiInfiniteDislocationElement(a2 + d2, x, x0, t, c);
  } else {
    CoefficientsForSemiInfiniteDislocationElement(a2 + d2, x, x1, t, c);
  }

  VectorProduct(d, b, dxb);
  VectorProduct(t, b, txb);
  TensorProduct(t, t, tct);
  TensorProduct(d, d, dcd);
  SymmetricTensorProduct(t, txb, tctxb_sym);
  SymmetricTensorProduct(t, d, tcd_sym);
  SymmetricTensorProduct(t, dxb, tcdxb_sym);
  SymmetricTensorProduct(d, txb, dctxb_sym);
  dxbdt = ScalarProduct(dxb, t);

  // Loop for symmetric
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = iDim; jDim < 3; jDim++) {
      stress[iDim][jDim] +=
          c[0] * (g4pv * dxbdt * I[iDim][jDim] + g4pv * dctxb_sym[iDim][jDim] -
                  g4p * tcdxb_sym[iDim][jDim]) +
          c[1] * (-1.0 * g4pv * v * tctxb_sym[iDim][jDim]) +
          c[2] * (g4pv * dxbdt * dcd[iDim][jDim] +
                  g4pv * dxbdt * a2 * I[iDim][jDim] -
                  g8p * a2 * tcdxb_sym[iDim][jDim]) +
          c[3] * (g8p * a2 * tctxb_sym[iDim][jDim] -
                  g4pv * dxbdt * tcd_sym[iDim][jDim]) +
          c[4] * (g4pv * dxbdt * tct[iDim][jDim]);
    }
  }
}
#endif

void DD_DislocationStress(double x[3], DD_t *dd, double stress[3][3]) {
  // Initialize
  InitializeTensor(stress);

  // Stress due to dislocation elements
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    FiniteDislocationElementStress(iElement, x, dd, stress);
  }

#ifdef _NOTERMINATION_
  // Stress due to semi-infinite dislocation element
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    for (int iNode = 0; iNode < 2; iNode++) {
      int nodeID = element->nodeID[iNode];
      NODE_t *node = &dd->node[nodeID];

      if (node->type != IMMOBILE_NODE && node->nElements == 1) {
        SemiInfiniteDislocationElementStress(iElement, x, dd, stress);
      }
    }
  }
#endif

  // Symmetric
  SymmetricTensor(stress);
}
