//========================================================================
// Reference:
//------------------------------------------------------------------------
// Authors:
// M. Lazar
// Title:
// A non-singular continuum theory of point defects using gradient 
// elasticity of bi-Helmholtz type
// Journal: Philosophical Magazine, Vol. 99, (2019), pp. 1563-1601
//-------------------------------------------------------------------------
// The implemented equation is for case (1) in the paper.
//=========================================================================
#include <math.h>
#include <stdio.h>

#include "incstr_math.h"
#include "incstr_struct.h"

static void DilatationCenterStress(double rx[3], double c, double a0,
                                   double stress[3][3]) {
  double c1 = 0.4 * a0;
  double c2 = 0.2 * a0;
  double r = INCSTR_VectorLength(3, rx);
  double r2 = r * r;
  double r3 = r2 * r;
  double r5 = r2 * r3;
  double e1 = exp(-1.0 * r / c1);
  double e2 = exp(-1.0 * r / c2);
  double f2 = 1.0 - 1.0 / (c1 * c1 - c2 * c2) * (c1 * c1 * e1 - c2 * c2 * e2) -
              r / (c1 * c1 - c2 * c2) * (c1 * e1 - c2 * e2) -
              r * r / (3.0 * (c1 * c1 - c2 * c2)) * (e1 - e2);
  double gbh = 1.0 / (4.0 * M_PI * (c1 * c1 - c2 * c2) * r) * (e1 - e2);
  double d[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  for (int iDim = 0; iDim < 3; iDim++) {
    double ri = rx[iDim];
    for (int jDim = 0; jDim < 3; jDim++) {
      double rj = rx[jDim];
      double dij = d[iDim][jDim];

      stress[iDim][jDim] = c * ((dij / r3 - 3.0 * ri * rj / r5) * f2 -
                                8.0 * M_PI / 3.0 * dij * gbh);
    }
  }
}

void INCSTR_DilatationCenterStress(double x[3], int id, INCLUSION_t *inclusion,
                                   MATERIAL_t *material, double stress[3][3]) {
  double *px = inclusion->x[id];
  double q = inclusion->eigen_strain[id];
  double a0 = inclusion->size[id][0];
  double g = material->g;
  double v = material->v;
  double c = g * q / (2.0 * M_PI) * (1.0 + v) / (1.0 - v);
  double **offset = material->offset;
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  for (int i = 0; i < material->nImages; i++) {
    double rx[3];
    for (int iDim = 0; iDim < 3; iDim++) {
      rx[iDim] = x[iDim] - (px[iDim] + offset[i][iDim]);
    }
    DilatationCenterStress(rx, c, a0, local_stress);
  }

  INCSTR_AddTensor(local_stress, stress);
}