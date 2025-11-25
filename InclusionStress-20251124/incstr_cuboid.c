//=======================================================================
// Reference:
//-----------------------------------------------------------------------
// Authors:
// Y.P. Chiu
// Title:
// On the Stress Field Due to Initial Strains in a Cuboid Surrounded by
// an Infinite Elastic Space
// Journal:
// Journal of Applied Mechanics (Transactions of ASME), Vol. 44,
// (1977), pp.587-590
//========================================================================
#include <math.h>
#include <stdio.h>

#include "incstr_math.h"
#include "incstr_struct.h"

static double Diiii(double ci, double cj, double ck, double r) {
  double t1 = atan((cj * ck) / (ci * r));
  double t2 = (ci * cj * ck) / (2.0 * r);
  double t3 = 1.0 / (ci * ci + cj * cj) + 1.0 / (ci * ci + ck * ck);

  return 2.0 * M_PI * M_PI * (t1 - t2 * t3);
}

static double Diiij(double ci, double cj, double ck, double r) {
  double sign = 1.0;
  double t1, t2;

  if (ck < 0.0) sign = -1.0;

  t1 = sign * log((r + fabs(ck)) / sqrt(ci * ci + cj * cj));
  t2 = (ci * ci * ck) / ((ci * ci + cj * cj) * r);

  return -1.0 * M_PI * M_PI * (t1 - t2);
}

static double Diijj(double ci, double cj, double ck, double r) {
  return (M_PI * M_PI * ci * cj * ck) / ((ci * ci + cj * cj) * r);
}

static double Diijk(double ci, double cj, double ck, double r) {
  return -1.0 * (M_PI * M_PI * ci) / r;
}

static void DTensor(double c[3], double d[3][3][3][3]) {
  double r = INCSTR_VectorLength(3, c);
  double c0 = c[0];
  double c1 = c[1];
  double c2 = c[2];

  d[0][0][0][0] = Diiii(c0, c1, c2, r);
  d[1][1][1][1] = Diiii(c1, c2, c0, r);
  d[2][2][2][2] = Diiii(c2, c0, c1, r);

  d[0][0][0][1] = d[0][0][1][0] = d[0][1][0][0] = d[1][0][0][0] =
      Diiij(c0, c1, c2, r);
  d[0][0][0][2] = d[0][0][2][0] = d[0][2][0][0] = d[2][0][0][0] =
      Diiij(c0, c2, c1, r);
  d[1][1][1][0] = d[1][1][0][1] = d[1][0][1][1] = d[0][1][1][1] =
      Diiij(c1, c0, c2, r);
  d[1][1][1][2] = d[1][1][2][1] = d[1][2][1][1] = d[2][1][1][1] =
      Diiij(c1, c2, c0, r);
  d[2][2][2][0] = d[2][2][0][2] = d[2][0][2][2] = d[0][2][2][2] =
      Diiij(c2, c0, c1, r);
  d[2][2][2][1] = d[2][2][1][2] = d[2][1][2][2] = d[1][2][2][2] =
      Diiij(c2, c1, c0, r);

  d[0][0][1][1] = d[0][1][0][1] = d[1][0][0][1] = d[0][1][1][0] =
      d[1][0][1][0] = d[1][1][0][0] = Diijj(c0, c1, c2, r);
  d[0][0][2][2] = d[0][2][0][2] = d[2][0][0][2] = d[0][2][2][0] =
      d[2][0][2][0] = d[2][2][0][0] = Diijj(c0, c2, c1, r);
  d[1][1][2][2] = d[1][2][1][2] = d[2][1][1][2] = d[1][2][2][1] =
      d[2][1][2][1] = d[2][2][1][1] = Diijj(c1, c2, c0, r);

  d[0][0][1][2] = d[0][1][0][2] = d[1][0][0][2] = d[0][1][2][0] =
      d[1][0][2][0] = d[1][2][0][0] = d[0][0][2][1] = d[0][2][0][1] =
          d[2][0][0][1] = d[0][2][1][0] = d[2][0][1][0] = d[2][1][0][0] =
              Diijk(c0, c1, c2, r);
  d[1][1][0][2] = d[1][0][1][2] = d[0][1][1][2] = d[1][0][2][1] =
      d[0][1][2][1] = d[0][2][1][1] = d[1][1][2][0] = d[1][2][1][0] =
          d[2][1][1][0] = d[1][2][0][1] = d[2][1][0][1] = d[2][0][1][1] =
              Diijk(c1, c0, c2, r);
  d[2][2][0][1] = d[2][0][2][1] = d[0][2][2][1] = d[2][0][1][2] =
      d[0][2][1][2] = d[0][1][2][2] = d[2][2][1][0] = d[2][1][2][0] =
          d[1][2][2][0] = d[2][1][0][2] = d[1][2][0][2] = d[1][0][2][2] =
              Diijk(c2, c0, c1, r);
}

static void DDisplacement(double g, double v, double lambda, double c[8][3],
                          double e[3][3], double ddisp[3][3]) {
  double d[8][3][3][3][3];

  for (int i = 0; i < 8; i++) {
    DTensor(c[i], d[i]);
  }

  for (int i = 0; i < 3; i++) {
    for (int q = 0; q < 3; q++) {
      double sign = 1.0;
      double v = 0.0;
      for (int n = 0; n < 8; n++) {
        double t1 = 1.0;
        double t2 = 1.0;
        double t3 = 1.0;
        double s1 = 0.0, s2 = 0.0;

        sign *= -1.0;

        t1 *= (1.0 - 2.0 * v) / (1.0 - v) * lambda;
        s1 = 0.0;
        for (int k = 0; k < 3; k++) {
          s1 += e[k][k];
        }
        t1 *= s1;
        s2 = 0.0;
        for (int m = 0; m < 3; m++) {
          s2 += d[n][i][q][m][m];
        }
        t1 *= s2;

        s2 = 0.0;
        t2 *= 4.0 * g;
        for (int j = 0; j < 3; j++) {
          s1 = 0.0;
          for (int m = 0; m < 3; m++) {
            s1 += d[n][j][q][m][m];
          }
          s2 += e[i][j] * s1;
        }
        t2 *= s2;

        s1 = 0.0;
        t3 *= (2.0 * g) / (1.0 - v);
        for (int j = 0; j < 3; j++) {
          for (int m = 0; m < 3; m++) {
            s1 += e[m][j] * d[n][i][q][m][j];
          }
        }
        t3 *= s1;

        v += sign * (t1 + t2 - t3);
      }
      ddisp[i][q] = v / (8.0 * M_PI * M_PI * M_PI * 2.0 * g);
    }
  }
}

static void CuboidalInclusionStress(double x[3], double size[3], double e,
                                    double g, double v,
                                    double rotation_tensor[3][3], int *inside,
                                    double stress[3][3]) {
  double rx[3];
  double lambda = (2.0 * v * g) / (1.0 - 2.0 * v);
  double y = 2.0 * g * (1.0 + v);
  double c = y / ((1.0 + v) * (1.0 - 2.0 * v));
  double c_tensor[8][3];
  double e_tensor[3][3] = {{e, 0.0, 0.0}, {0.0, e, 0.0}, {0.0, 0.0, e}};
  double ddisp[3][3];
  double s_tensor[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  INCSTR_RotateCoordinateVector(rotation_tensor, x, rx);

  c_tensor[0][0] = rx[0] - size[0];
  c_tensor[0][1] = rx[1] - size[1];
  c_tensor[0][2] = rx[2] - size[2];
  c_tensor[1][0] = rx[0] + size[0];
  c_tensor[1][1] = rx[1] - size[1];
  c_tensor[1][2] = rx[2] - size[2];
  c_tensor[2][0] = rx[0] + size[0];
  c_tensor[2][1] = rx[1] + size[1];
  c_tensor[2][2] = rx[2] - size[2];
  c_tensor[3][0] = rx[0] - size[0];
  c_tensor[3][1] = rx[1] + size[1];
  c_tensor[3][2] = rx[2] - size[2];
  c_tensor[4][0] = rx[0] - size[0];
  c_tensor[4][1] = rx[1] + size[1];
  c_tensor[4][2] = rx[2] + size[2];
  c_tensor[5][0] = rx[0] - size[0];
  c_tensor[5][1] = rx[1] - size[1];
  c_tensor[5][2] = rx[2] + size[2];
  c_tensor[6][0] = rx[0] + size[0];
  c_tensor[6][1] = rx[1] - size[1];
  c_tensor[6][2] = rx[2] + size[2];
  c_tensor[7][0] = rx[0] + size[0];
  c_tensor[7][1] = rx[1] + size[1];
  c_tensor[7][2] = rx[2] + size[2];

  DDisplacement(g, v, lambda, c_tensor, e_tensor, ddisp);
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      if (fabs(rx[0]) < size[0] && fabs(rx[1]) < size[1] &&
          fabs(rx[2]) < size[2]) {
        s_tensor[iDim][jDim] = 0.5 * (ddisp[iDim][jDim] + ddisp[jDim][iDim]) -
                               e_tensor[iDim][jDim];
        (*inside) = 1;
      } else {
        s_tensor[iDim][jDim] = 0.5 * (ddisp[iDim][jDim] + ddisp[jDim][iDim]);
      }
    }
  }

  local_stress[0][0] =
      c * ((1.0 - v) * s_tensor[0][0] + v * (s_tensor[1][1] + s_tensor[2][2]));
  local_stress[1][1] =
      c * ((1.0 - v) * s_tensor[1][1] + v * (s_tensor[2][2] + s_tensor[0][0]));
  local_stress[2][2] =
      c * ((1.0 - v) * s_tensor[2][2] + v * (s_tensor[0][0] + s_tensor[1][1]));
  local_stress[0][1] = local_stress[1][0] = 2.0 * g * s_tensor[0][1];
  local_stress[1][2] = local_stress[2][1] = 2.0 * g * s_tensor[1][2];
  local_stress[2][0] = local_stress[0][2] = 2.0 * g * s_tensor[2][0];

  INCSTR_RotateStressTensor(rotation_tensor, local_stress, stress);
}

static void RotationTensor(double **matrix_tensor, double **inclusion_tensor,
                           double tensor[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double v = 0.0;
      for (int k = 0; k < 3; k++) {
        v += matrix_tensor[i][k] * inclusion_tensor[j][k];
      }
      tensor[i][j] = v;
    }
  }
}

void INCSTR_CuboidalInclusionStress(double x[3], int id, INCLUSION_t *inclusion,
                                    MATERIAL_t *material, int *inside,
                                    double stress[3][3]) {
  double e = inclusion->eigen_strain[id];
  double *px = inclusion->x[id];
  double size[3] = {inclusion->size[id][0], inclusion->size[id][1],
                    inclusion->size[id][2]};
  double g = material->g;
  double v = material->v;
  double **offset = material->offset;
  double rotation_tensor[3][3];
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  RotationTensor(material->orientation, inclusion->orientation[id],
                 rotation_tensor);
  for (int i = 0; i < material->nImages; i++) {
    double rx[3];

    for (int iDim = 0; iDim < 3; iDim++) {
      rx[iDim] = x[iDim] - (px[iDim] + offset[i][iDim]);
    }
    CuboidalInclusionStress(rx, size, e, g, v, rotation_tensor, inside,
                            local_stress);
  }

  INCSTR_AddTensor(local_stress, stress);
}