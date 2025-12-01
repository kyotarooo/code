//========================================================================
// Reference:
//------------------------------------------------------------------------
// Authors:
// A.L. Kolenikova, M.Yu. Gutkin, A.E. Romanov
// Title:
// Analytical elastic models of finite cylindrical and truncated spherical
// inclusions
// Journal:
// International Journal of Solids and Structures, Vol. 143,
// (2018), pp.59-72
//=========================================================================
#include <math.h>
#include <stdio.h>

#include "incstr_bessel.h"
#include "incstr_elliptic.h"
#include "incstr_math.h"
#include "incstr_struct.h"

static void RotationTensor(double rx[3], double ma[3], double tensor[3][3]) {
  double d = 0.0;
  double r[3], p[3], z[3];

  for (int i = 0; i < 3; i++) {
    z[i] = ma[i];
    d += rx[i] * z[i];
  }
  for (int i = 0; i < 3; i++) {
    r[i] = rx[i] - d * z[i];
  }
  if (d < 0.0) {
    for (int i = 0; i < 3; i++) {
      z[i] *= -1.0;
    }
  }

  d = INCSTR_VectorLength(3, r);
  if (d < 1.0e-32) {
    double vr[3] = {0.413412124, -0.2414214, 0.2314211};  // Random direction

    INCSTR_VectorProduct(vr, z, r);
    d = INCSTR_VectorLength(3, r);
  }
  for (int i = 0; i < 3; i++) {
    r[i] /= d;
  }

  INCSTR_VectorProduct(z, r, p);

  for (int i = 0; i < 3; i++) {
    tensor[0][i] = r[i];
    tensor[1][i] = p[i];
    tensor[2][i] = z[i];
  }
}

static void MoveGridPointOutsideBoundaryLayer(double rx[3], double c, double b,
                                              double thickness) {
  double rd, rz;
  double dd = 0.0, dz = 0.0;
  double eps = thickness * 1.0e-06;

  rd = INCSTR_VectorLength(2, rx);
  rz = rx[2];

  if ((fabs(rd - c) < thickness) && (fabs(rz) < b + thickness)) {
    if (rd > c) {
      dd += thickness - (rd - c);
    } else {
      dd -= thickness - (c - rd);
    }
  }
  if ((fabs(fabs(rz) - b) < thickness) && (fabs(rd) < c + thickness)) {
    if (fabs(rz) > b) {
      dz += thickness - (fabs(rz) - b);
    } else {
      dz -= thickness - (b - fabs(rz));
    }
  }
  if (dd > 0.0 && dz > 0.0) {
    if (dd < dz)
      dz = 0.0;
    else
      dd = 0.0;
  }

  if (rd > eps) {
    for (int iDim = 0; iDim < 2; iDim++) {
      rx[iDim] += dd * rx[iDim] / rd;
    }
  }
  if (fabs(rx[2]) > eps) {
    rx[2] += dz * rx[2] / fabs(rx[2]);
  }
}

static void CylindricalInclusionStress(double x[3], double r, double h,
                                       double thickness, double ma[3], double c,
                                       int* inside, double stress[3][3]) {
  double rotation_tensor[3][3];
  double rx[3];
  double rd = 0.0, rz;
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double r_, sqr_, h_, m1, m2, m3, x1, x2, x3;
  double srr, spp, szz, srz;
  double (*BF)(int index1, int index2, int index3, double param[7]) =
      INCSTR_BesselFunction;
  double (*Ek)(double m) = INCSTR_EllipticIntegralEk;
  double (*Fk)(double m) = INCSTR_EllipticIntegralFk;
  double (*Pik)(double h, double m) = INCSTR_EllipticIntegralPik;

  RotationTensor(x, ma, rotation_tensor);
  INCSTR_RotateCoordinateVector(rotation_tensor, x, rx);
  MoveGridPointOutsideBoundaryLayer(rx, r, 0.5 * h, thickness);
  rd = INCSTR_VectorLength(2, rx);
  rz = rx[2];

  r_ = rd / r;
  sqr_ = sqrt(r_);
  h_ = 4.0 * r_ / ((1.0 + r_) * (1.0 + r_));

  x1 = (0.5 * h - fabs(rx[2])) / r;
  x2 = (0.5 * h + fabs(rx[2])) / r;
  x3 = (-0.5 * h + fabs(rx[2])) / r;
  m1 = sqrt(4.0 * r_ / ((1.0 + r_) * (1.0 + r_) + x1 * x1));
  m2 = sqrt(4.0 * r_ / ((1.0 + r_) * (1.0 + r_) + x2 * x2));
  m3 = sqrt(4.0 * r_ / ((1.0 + r_) * (1.0 + r_) + x3 * x3));

  if (rz >= -0.5 * h && rz <= 0.5 * h) {
    double param1[7] = {r_, sqr_, x1, m1, Ek(m1), Fk(m1), Pik(h_, m1)};
    double param2[7] = {r_, sqr_, x2, m2, Ek(m2), Fk(m2), Pik(h_, m2)};

    if (rd <= r) {
      srr = -0.5 * c *
            (BF(1, 0, 0, param1) + BF(1, 0, 0, param2) - BF(1, 2, 0, param1) -
             BF(1, 2, 0, param2) + 2.0);
      spp = -1.0 * c *
            (r / rd * (BF(1, 1, -1, param1) + BF(1, 1, -1, param2)) + 1.0);
      szz = c * (BF(1, 0, 0, param1) + BF(1, 0, 0, param2) - 2.0);
      if (rz > 0.0) {
        srz = -1.0 * c * (BF(1, 1, 0, param1) - BF(1, 1, 0, param2));
      } else {
        srz = -1.0 * c * (-1.0) * (BF(1, 1, 0, param1) - BF(1, 1, 0, param2));
      }
      (*inside) = 1;
    } else {
      srr = 0.5 * c *
            ((-2.0 * r * r) / (rd * rd) - BF(1, 0, 0, param1) -
             BF(1, 0, 0, param2) + BF(1, 2, 0, param1) + BF(1, 2, 0, param2));
      spp =
          c * (r / rd * (r / rd - BF(1, 1, -1, param1) - BF(1, 1, -1, param2)));
      szz = c * (BF(1, 0, 0, param1) + BF(1, 0, 0, param2));
      if (rz > 0.0) {
        srz = -1.0 * c * (BF(1, 1, 0, param1) - BF(1, 1, 0, param2));
      } else {
        srz = -1.0 * c * (-1.0) * (BF(1, 1, 0, param1) - BF(1, 1, 0, param2));
      }
    }
  } else {
    double param2[7] = {r_, sqr_, x2, m2, Ek(m2), Fk(m2), Pik(h_, m2)};
    double param3[7] = {r_, sqr_, x3, m3, Ek(m3), Fk(m3), Pik(h_, m3)};

    srr = 0.5 * c *
          (BF(1, 0, 0, param3) - BF(1, 0, 0, param2) - BF(1, 2, 0, param3) +
           BF(1, 2, 0, param2));
    spp = c * (r / rd * (BF(1, 1, -1, param3) - BF(1, 1, -1, param2)));
    szz = c * (-1.0 * BF(1, 0, 0, param3) + BF(1, 0, 0, param2));
    if (rz > 0.0) {
      srz = -1.0 * c * (BF(1, 1, 0, param3) - BF(1, 1, 0, param2));
    } else {
      srz = -1.0 * c * (-1.0) * (BF(1, 1, 0, param3) - BF(1, 1, 0, param2));
    }
  }

  local_stress[0][0] = srr;
  local_stress[1][1] = spp;
  local_stress[2][2] = szz;
  local_stress[0][2] = local_stress[2][0] = srz;

  INCSTR_RotateStressTensor(rotation_tensor, local_stress, stress);
}

void INCSTR_CylindricalInclusionStress(double x[3], int id,
                                       INCLUSION_t* inclusion,
                                       MATERIAL_t* material, int* inside,
                                       double stress[3][3]) {
  double e = inclusion->eigen_strain[id];
  double* px = inclusion->x[id];
  double r = inclusion->size[id][0];
  double h = inclusion->size[id][1];
  double thickness = 0.5 * inclusion->size[id][2];
  double g = material->g;
  double v = material->v;
  double c = g * (1.0 + v) * e / (1.0 - v);
  double ma[3];
  double* major_axis = inclusion->orientation[id][0];
  double** offset = material->offset;
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  for (int iDim = 0; iDim < 3; iDim++) {
    double v = 0.0;
    for (int jDim = 0; jDim < 3; jDim++) {
      v += material->orientation[iDim][jDim] * major_axis[jDim];
    }
    ma[iDim] = v;
  }

  for (int i = 0; i < material->nImages; i++) {
    double rx[3];
    for (int iDim = 0; iDim < 3; iDim++) {
      rx[iDim] = x[iDim] - (px[iDim] + offset[i][iDim]);
    }
    CylindricalInclusionStress(rx, r, h, thickness, ma, c, inside,
                               local_stress);
  }

  INCSTR_AddTensor(local_stress, stress);
}