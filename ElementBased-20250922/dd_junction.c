#include <math.h>
#include <stdio.h>

#include "dd_defs.h"
#include "dd_element_handle.h"
#include "dd_periodic_bc.h"
#include "dd_struct.h"

static void SubstituteVector(double *a, double b[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] = a[iDim];
  }
}

static void VectorProduct(double a[3], double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static double ScalarProduct(double a[3], double b[3]) {
  double c = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    c += a[iDim] * b[iDim];
  }
  return c;
}

static void UnitVector(double a[3]) {
  double b = sqrt(ScalarProduct(a, a));

  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] /= b;
  }
}

static void UpdatePosition(double a, double b[3], double x[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    x[iDim] += a * b[iDim];
  }
}

void DD_Junction(int elementID0, int elementID1, DD_t *dd) {
  ELEMENT_t *element0 = &dd->element[elementID0];
  ELEMENT_t *element1 = &dd->element[elementID1];
  NODE_t *node00 = &dd->node[element0->nodeID[0]];
  NODE_t *node01 = &dd->node[element0->nodeID[1]];
  NODE_t *node10 = &dd->node[element1->nodeID[0]];
  NODE_t *node11 = &dd->node[element1->nodeID[1]];
  double x00[3], x01[3], x10[3], x11[3];
  double b0[3], b1[3];
  double t0[3], t1[3], t[3];
  double s0[3], s1[3];
  double g0[3], g1[3];
  double d0, d1;
  double sg01, sg10;
  double d;

  SubstituteVector(node00->x, x00);
  SubstituteVector(node01->x, x01);
  SubstituteVector(node10->x, x10);
  SubstituteVector(node11->x, x11);
  DD_PeriodicReturn(x00, x01, &dd->bc);
  DD_PeriodicReturn(x00, x10, &dd->bc);
  DD_PeriodicReturn(x10, x11, &dd->bc);

  SubstituteVector(element0->burgers, b0);
  SubstituteVector(element1->burgers, b1);
  SubstituteVector(element0->tangent, t0);
  SubstituteVector(element1->tangent, t1);
  SubstituteVector(element0->slip, s0);
  SubstituteVector(element1->slip, s1);

  VectorProduct(s0, t0, g0);
  VectorProduct(s1, t1, g1);

  d0 = ScalarProduct(s0, x00);
  d1 = ScalarProduct(s1, x10);
  sg01 = ScalarProduct(s0, g1);
  sg10 = ScalarProduct(s1, g0);

  d = (d1 - ScalarProduct(s1, x00)) / sg10;
  UpdatePosition(d, g0, x00);
  d = (d1 - ScalarProduct(s1, x01)) / sg10;
  UpdatePosition(d, g0, x01);
  d = (d0 - ScalarProduct(s0, x10)) / sg01;
  UpdatePosition(d, g1, x10);
  d = (d0 - ScalarProduct(s0, x11)) / sg01;
  UpdatePosition(d, g1, x11);

  // Merge elements to make a stair-rod dislocation element of junction
  if (ScalarProduct(t0, t1) < 0.0) {
    for (int iDim = 0; iDim < 3; iDim++) {
      node00->x[iDim] = 0.5 * (x00[iDim] + x11[iDim]);
      node01->x[iDim] = 0.5 * (x01[iDim] + x10[iDim]);
      element0->burgers[iDim] -= b1[iDim];
    }
  } else {
    for (int iDim = 0; iDim < 3; iDim++) {
      node00->x[iDim] = 0.5 * (x00[iDim] + x10[iDim]);
      node01->x[iDim] = 0.5 * (x01[iDim] + x11[iDim]);
      element0->burgers[iDim] += b1[iDim];
    }
  }

  VectorProduct(s0, s1, t);
  UnitVector(t);
  node00->type = CONSTRAINT_NODE;
  node01->type = CONSTRAINT_NODE;
  for (int iDim = 0; iDim < 3; iDim++) {
    node00->direction[iDim] = t[iDim];
    node01->direction[iDim] = t[iDim];
    element0->slip[iDim] = t[iDim];
  }
  DD_MergeElements(elementID0, elementID1, dd);
}