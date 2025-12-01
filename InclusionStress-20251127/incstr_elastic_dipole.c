//==============================================================================
// Reference:
//------------------------------------------------------------------------------
// Authors:
// G. Gengor, O.K. Celebi, A.S.K. Mohammed, H. Sehitoglu
// Title:
// Continuum strain of point defects
// Journal: Journal of the Mechanics and Physics of Solids, Vol. 188, (2024), p.
// 105653
//==============================================================================
// To make it non-singular
//------------------------------------------------------------------------------
// Authors:
// W. Cai, A. Arsenlis, C.R. Weinberger, V.V. Bulatov
// Title:
// A non-singular continuum theory of dislocations
// Journal: Journal of the Mechanics and Physics of Solids, Vol. 54, (2006), pp.
// 561-587
//===============================================================================

#include <math.h>
#include <stdio.h>

#include "incstr_math.h"
#include "incstr_struct.h"

static void SecondDerivativeGreenFunction(double x[3], double a,
                                          MATERIAL_t* material,
                                          double ddg[3][3][3][3]) {
  double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + a * a);
  double g = material->g;
  double v = material->v;
  double c0 = 1.0 / (16.0 * M_PI * g * (1.0 - v) * r * r * r);
  double c1 = 3.0 - 4.0 * v;
  double d[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  double rx[3];

  for (int i = 0; i < 3; i++) {
    rx[i] = x[i] / r;
  }

  for (int i = 0; i < 3; i++) {
    double ri = rx[i];
    for (int j = 0; j < 3; j++) {
      double dij = d[i][j];
      double rj = rx[j];
      double rij = ri * rj;
      for (int k = 0; k < 3; k++) {
        double dik = d[i][k];
        double djk = d[j][k];
        double rk = rx[k];
        double rik = ri * rk;
        double rjk = rj * rk;
        for (int l = 0; l < 3; l++) {
          double dil = d[i][l];
          double djl = d[j][l];
          double dkl = d[k][l];
          double rl = rx[l];
          double ril = ri * rl;
          double rjl = rj * rl;
          double rkl = rk * rl;

          ddg[i][j][k][l] =
              c0 * (dik * djl + dil * djk - c1 * dij * (dkl - 3.0 * rkl) -
                    3.0 * (dkl * rij + dik * rjl + djk * ril + dil * rjk +
                           djl * rik) +
                    15.0 * rij * rkl);
        }
      }
    }
  }
}

static void ElasticDipoleStress(double x[3], double c[3][3][3][3], double a,
                                MATERIAL_t* material,
                                double elastic_dipole[3][3],
                                double stress[3][3]) {
  double ddg[3][3][3][3];
  double du[3][3];

  SecondDerivativeGreenFunction(x, a, material, ddg);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double v = 0.0;
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          v += ddg[i][k][l][j] * elastic_dipole[k][l];
        }
      }
      du[i][j] = v;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double v = 0.0;
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          v += c[i][j][k][l] * du[k][l];
        }
      }
      stress[i][j] += v;
    }
  }
}

void INCSTR_ElasticDipoleStress(double x[3], int id, INCLUSION_t* inclusion,
                                MATERIAL_t* material, double stress[3][3]) {
  double* px = inclusion->x[id];
  double a = inclusion->size[id][0];
  double elastic_dipole[3][3];
  int iidx[9] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
  int jidx[9] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  double** offset = material->offset;
  double g = material->g;
  double v = material->v;
  double lambda = 2.0 * g * v / (1.0 - 2.0 * v);
  double c[3][3][3][3];
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  for (int i = 0; i < 9; i++) {
    elastic_dipole[iidx[i]][jidx[i]] = inclusion->elastic_dipole[id][i];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          c[i][j][k][l] = 0.0;
        }
      }
    }
  }
  c[0][0][0][0] = c[1][1][1][1] = c[2][2][2][2] = lambda + 2.0 * g;
  c[0][0][1][1] = c[0][0][2][2] = c[1][1][0][0] = c[1][1][2][2] =
      c[2][2][0][0] = c[2][2][1][1] = lambda;
  c[0][1][0][1] = c[0][1][1][0] = c[1][0][0][1] = c[1][0][1][0] =
      c[1][2][1][2] = c[1][2][2][1] = c[2][1][1][2] = c[2][1][2][1] =
          c[0][2][0][2] = c[0][2][2][0] = c[2][0][0][2] = c[2][0][2][0] = g;

  for (int i = 0; i < material->nImages; i++) {
    double rx[3];
    for (int iDim = 0; iDim < 3; iDim++) {
      rx[iDim] = x[iDim] - (px[iDim] + offset[i][iDim]);
    }
    ElasticDipoleStress(rx, c, a, material, elastic_dipole, local_stress);
  }

  INCSTR_AddTensor(local_stress, stress);
}