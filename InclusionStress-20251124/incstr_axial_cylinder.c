//========================================================================
// Reference:
//------------------------------------------------------------------------
// Authors:
// H. Hasegawa, V.G. Lee, T. Mura
// Title:
// The stress fields caused by a circular cylindrical inclusion
// Journal:
// Journal of Applied Mechanics, Vol. 59, (1992), pp.107-114
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

static void AxialCylindricalInclusionStress(double x[3], double c, double b,
                                            double thickness, double ma[3],
                                            double e, MATERIAL_t* material,
                                            int* inside, double stress[3][3]) {
  double rotation_tensor[3][3];
  double rx[3], rd, rz, rd2, rdoc, crd, cprd, cprd2, cmrd, cmrd2, cmrdocprd;
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double srr = 0.0, spp = 0.0, szz = 0.0, srz = 0.0;
  double g = material->g;
  double v = material->v;
  double cc12 = 1.0 - 2.0 * v;
  double cc1 = 1.0 - v;
  double c2 = c * c;
  double sign[2] = {-1.0, 1.0};

  RotationTensor(x, ma, rotation_tensor);
  INCSTR_RotateCoordinateVector(rotation_tensor, x, rx);
  MoveGridPointOutsideBoundaryLayer(rx, c, b, thickness);
  rd = INCSTR_VectorLength(2, rx);

  rd2 = rd * rd;
  rdoc = rd / c;
  crd = c * rd;
  cprd = c + rd;
  cprd2 = cprd * cprd;
  cmrd = c - rd;
  cmrd2 = cmrd * cmrd;
  cmrdocprd = cmrd / (c + rd);
  rz = rx[2];
  if (fabs(rd) < c && fabs(rz) < b) {
    (*inside) = 1;
  }

  double z[2] = {rz - b, rz + b};
  double cc = g * c * e / (4.0 * M_PI * sqrt(crd) * cc1 * cc12);
  double p = 4.0 * crd / cprd2;

  for (int i = 0; i < 2; i++) {
    double sbbrr1, sbbrr3, sbbpp1, sbbpp3, sbbzz1, sbbzz3, sbbrz1, sbbrz3;
    double zi = z[i];
    double zi2 = zi * zi;
    double k = sqrt(4.0 * crd / (cprd2 + zi2));
    double xi = (rd2 + c2 + zi2) / (2.0 * crd);
    double xi2 = xi * xi;
    double xi21inv = 1.0 / (xi2 - 1.0);
    double sqrtx1 = sqrt(xi + 1.0);
    double kx = M_SQRT2 / sqrtx1;
    double fi = INCSTR_EllipticIntegralFk(kx);
    double ei = INCSTR_EllipticIntegralEk(kx);
    double kpik = k * INCSTR_EllipticIntegralPik(p, k);
    double qn = kx * fi;
    double qp = xi * kx * fi - M_SQRT2 * sqrtx1 * ei;
    double gi = xi * qp - qn;
    double xi21invgi = xi21inv * gi;
    double signi = sign[i];

    sbbrr1 = -1.0 * zi / rd *
             ((3.0 - 4.0 * v) * qp - (rd2 - cc12 * c2) / crd * qn +
              (xi - rdoc) * xi21invgi -
              (2.0 * cc1 * rdoc * cmrdocprd + cc12 * cmrd2 / crd) * kpik);
    sbbpp1 = zi / rd *
             (4.0 * cc1 * qp + (rd2 + cc12 * c2) / crd * qn +
              (2.0 * v * rdoc * cmrdocprd - cc12 * cmrd2 / crd) * kpik);
    sbbzz1 = zi / rd *
             (2.0 * v * rdoc * (qn + cmrdocprd * kpik) - qp +
              (xi - rdoc) * xi21invgi);
    sbbrz1 = 2.0 * cc12 * qp + zi2 / crd * xi21invgi;
    sbbrr3 = -1.0 * zi / c *
             (2.0 * v * (qn + cmrdocprd * kpik) -
              c / rd * (qp + (xi - rdoc) * xi21invgi));
    sbbpp3 = -2.0 * zi / c * (c / rd * qp + v * (qn + cmrdocprd * kpik));
    sbbzz3 = -1.0 * zi / c *
             (2.0 * cc1 * (qn + cmrdocprd * kpik) -
              c / rd * (qp - (xi - rdoc) * xi21invgi));
    sbbrz3 = 2.0 * cc12 * qp - zi2 / crd * xi21invgi;

    srr += signi * (v * sbbrr1 + cc1 * sbbrr3);
    spp += signi * (v * sbbpp1 + cc1 * sbbpp3);
    szz += signi * (v * sbbzz1 + cc1 * sbbzz3);
    srz += signi * (v * sbbrz1 + cc1 * sbbrz3);
  }
  srr *= cc;
  spp *= cc;
  szz *= cc;
  srz *= cc;

  local_stress[0][0] = srr;
  local_stress[1][1] = spp;
  local_stress[2][2] = szz;
  local_stress[0][2] = local_stress[2][0] = srz;
  INCSTR_RotateStressTensor(rotation_tensor, local_stress, stress);
}

void INCSTR_AxialCylindricalInclusionStress(double x[3], int id,
                                            INCLUSION_t* inclusion,
                                            MATERIAL_t* material, int* inside,
                                            double stress[3][3]) {
  double e = inclusion->eigen_strain[id];
  double* px = inclusion->x[id];
  double c = inclusion->size[id][0];
  double b = 0.5 * inclusion->size[id][1];
  double thickness = 0.5 * inclusion->size[id][2];
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

    AxialCylindricalInclusionStress(rx, c, b, thickness, ma, e, material,
                                    inside, local_stress);
  }

  INCSTR_AddTensor(local_stress, stress);
}
