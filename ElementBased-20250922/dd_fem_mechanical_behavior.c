#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_periodic_bc.h"
#include "dd_struct.h"
#include "fem_struct.h"

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

static double VectorLength(double a[3]) {
  double r = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    r += a[iDim] * a[iDim];
  }
  return sqrt(r);
}

static double ScalarProduct(double a[3], double b[3]) {
  double c = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    c += a[iDim] * b[iDim];
  }
  return c;
}

static double SweptArea(double x00[3], double x01[3], double x10[3],
                        double x11[3], double n[3]) {
  double d0[3], d1[3];

  // Diagonal vector
  for (int iDim = 0; iDim < 3; iDim++) {
    d0[iDim] = x11[iDim] - x00[iDim];
    d1[iDim] = x01[iDim] - x10[iDim];
  }

  // Swept area is calculated with the vector product
  VectorProduct(d0, d1, n);
  return 0.5 * VectorLength(n);
}

static double SimulationVolume(double s[3], FEM_t *fem) {
  double volume = 1.0;

  if (fem->nElements == 0) {
    for (int iDim = 0; iDim < 3; iDim++) {
      volume *= s[iDim];
    }
  } else {
    volume = fem->volume;
  }
  return volume;
}

static void ElementPlasticStrainIncrement(int elementID, double volume,
                                          DD_t *dd,
                                          double dplastic_strain[3][3]) {
  ELEMENT_t *element = &dd->element[elementID];
  int nodeID0 = element->nodeID[0];
  int nodeID1 = element->nodeID[1];
  NODE_t *node0 = &dd->node[nodeID0];
  NODE_t *node1 = &dd->node[nodeID1];
  double b[3], s[3];
  double n[3];
  double x00[3], x01[3], x10[3], x11[3];
  double area;
  double sdn;

  // Position vectors before and after the time integral
  SubstituteVector(node0->x0, x00);
  SubstituteVector(node0->x, x01);
  SubstituteVector(node1->x0, x10);
  SubstituteVector(node1->x, x11);
  DD_PeriodicReturn(x00, x01, &dd->bc);
  DD_PeriodicReturn(x00, x10, &dd->bc);
  DD_PeriodicReturn(x00, x11, &dd->bc);

  // Area swept by the dislocation element
  area = SweptArea(x00, x01, x10, x11, n);

  // Element information
  SubstituteVector(element->burgers, b);
  SubstituteVector(element->slip, s);

  // Check glide direction
  sdn = ScalarProduct(s, n);
  if (sdn < 0.0) area *= -1.0;

  // Plastic strain increment
  area /= volume;
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      dplastic_strain[iDim][jDim] +=
          0.5 * area * (b[iDim] * s[jDim] + b[jDim] * s[iDim]);
    }
  }
}

void DD_FEM_MechanicalBehavior(DD_t *dd, FEM_t *fem) {
  double dstrain[3][3];
  double dplastic_strain[3][3];
  double delastic_strain[3][3];
  double dt = dd->step.dt;
  double volume = SimulationVolume(dd->bc.size, fem);
  int i = dd->mechanical_behavior.index[0];
  int j = dd->mechanical_behavior.index[1];
  double modulus;

  // Initialize
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      dstrain[iDim][jDim] =
          dd->mechanical_behavior.strain_rate[iDim][jDim] * dt;
      dplastic_strain[iDim][jDim] = 0.0;
    }
  }

  // Plastic strain increement of all elements
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ElementPlasticStrainIncrement(iElement, volume, dd, dplastic_strain);
  }

  // Update strains
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      delastic_strain[iDim][jDim] =
          dstrain[iDim][jDim] - dplastic_strain[iDim][jDim];
      dd->mechanical_behavior.strain[iDim][jDim] += dstrain[iDim][jDim];
      dd->mechanical_behavior.plastic_strain[iDim][jDim] +=
          dplastic_strain[iDim][jDim];
    }
  }

  // No stress update
  if (i == -1 || j == -1) return;

  // Elastic constant
  if (i == j)
    modulus = dd->material.e;
  else
    modulus = 2.0 * dd->material.g;

  // Update applied stress
  dd->bc.stress[i][j] += modulus * (dstrain[i][j] - dplastic_strain[i][j]);
  dd->bc.stress[j][i] = dd->bc.stress[i][j];

  // Update applied stress to FEM
  // fem->equation.scale_f = dd->bc.stress[i][j];
}
