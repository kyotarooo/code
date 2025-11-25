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
#include <stdlib.h>

#include "incstr_legendre.h"
#include "incstr_math.h"
#include "incstr_struct.h"

static void RotationTensor(double rx[3], double ma[3], double tensor[3][3]) {
  double d = 0.0;
  double r[3], t[3], p[3];

  for (int i = 0; i < 3; i++) {
    r[i] = rx[i];
  }
  d = INCSTR_VectorLength(3, r);
  if (d < 1.0e-20) {
    double vr[3] = {0.413412124, -0.2414214, 0.2314211};  // Random direction

    for (int i = 0; i < 3; i++) {
      r[i] = vr[i];
    }
    d = INCSTR_VectorLength(3, r);
  }
  for (int i = 0; i < 3; i++) {
    r[i] /= d;
  }

  INCSTR_VectorProduct(r, ma, p);
  d = INCSTR_VectorLength(3, p);
  if (d < 1.0e-20) {
    double vr[3] = {-0.113412124, -0.26614, 0.8314211};  // Random direction

    INCSTR_VectorProduct(r, vr, p);
    d = INCSTR_VectorLength(3, p);
  }
  for (int i = 0; i < 3; i++) {
    p[i] /= d;
  }

  INCSTR_VectorProduct(r, p, t);
  for (int i = 0; i < 3; i++) {
    tensor[0][i] = r[i];
    tensor[1][i] = t[i];
    tensor[2][i] = p[i];
  }
}

static void TruncatedSphericalInclusionStress(double x[3], double r, double z1,
                                              double z2, double ma[3], double c,
                                              int *inside,
                                              double stress[3][3]) {
  double rotation_tensor[3][3];
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double srr, spp, stt, srt;
  double z = 0.0, d = 0.0;
  double cs0, s0, ct0;
  double cs1 = z1 / r;
  double cs2 = z2 / r;
  double rr0;
  double p_0[NMAX_LEGENDRE + 1], p_1[NMAX_LEGENDRE + 1], p_2[NMAX_LEGENDRE + 1];
  double p1[NMAX_LEGENDRE + 1];

  RotationTensor(x, ma, rotation_tensor);
  for (int iDim = 0; iDim < 3; iDim++) {
    z += x[iDim] * ma[iDim];
  }
  d = INCSTR_VectorLength(3, x);
  rr0 = d / r;

  if (d < 1.0e-12 * r) {
    cs0 = 1.0;
  } else {
    cs0 = z / d;
  }
  s0 = sqrt(1.0 - cs0*cs0);
  if (s0 < 1.0e-12) {
    ct0 = 0.0;
  } else {
    ct0 = cs0 / s0;
  }
  
  // Legendre polynomials are calculated in advance.
  // The results are stored in arrays.
  INCSTR_LegendrePolynomials(cs0, p_0);
  INCSTR_LegendrePolynomials(cs1, p_1);
  INCSTR_LegendrePolynomials(cs2, p_2);
  INCSTR_AssociatedLegendrePolynomials(cs0, p1);

  if (d < r) {
    double rr0n = 1.0;

    srr = srt = stt = spp = 0.0;
    for (int i = 2; i <= NMAX_LEGENDRE; i++) {
      double bn = INCSTR_LureBn(i, p_1, p_2);
      double pn = p_0[i];
      double p1n = p1[i];
      double n = (double)i;

      srr += bn * n * (n - 1.0) * rr0n * pn;
      srt += bn * (n - 1.0) * rr0n * p1n;
      stt += bn * rr0n * (-1.0 * n * n * pn - ct0 * p1n);
      spp += bn * rr0n * (n * pn + ct0 * p1n);
      rr0n *= rr0;
    }
    if (z >= z1 && z <= z2) {
      srr += -2.0 + 2.0 * cs0 * cs0;
      srt += -2.0 * s0 * cs0;
      stt += -2.0 + 2.0 * s0 * s0;
      spp += -2.0;
      (*inside) = 1;
    }
  } else {
    double rr0n = 1.0 / (rr0 * rr0 * rr0);

    srr = srt = stt = spp = 0.0;
    for (int i = 0; i <= NMAX_LEGENDRE - 3; i++) {
      double dn = INCSTR_LureDn(i, p_1, p_2);
      double pn = p_0[i];
      double p1n = p1[i];
      double n = (double)i;

      srr += dn * (n + 1.0) * (n + 2.0) * rr0n * pn;
      srt += -1.0 * dn * (n + 2.0) * rr0n * p1n;
      stt += dn * rr0n * (-1.0 * (n + 1.0) * (n + 1.0) * pn - ct0 * p1n);
      spp += dn * rr0n * (-1.0 * (n + 1.0) * pn + ct0 * p1n);
      rr0n /= rr0;
    }
  }

  local_stress[0][0] = c * srr;
  local_stress[1][1] = c * stt;
  local_stress[2][2] = c * spp;
  local_stress[0][1] = local_stress[1][0] = c * srt;

  INCSTR_RotateStressTensor(rotation_tensor, local_stress, stress);
}

void INCSTR_TruncatedSphericalInclusionStress(double x[3], int id,
                                              INCLUSION_t *inclusion,
                                              MATERIAL_t *material, int *inside,
                                              double stress[3][3]) {
  double e = inclusion->eigen_strain[id];
  double *px = inclusion->x[id];
  double r = inclusion->size[id][0];
  double z1 = inclusion->size[id][1];
  double z2 = inclusion->size[id][2];
  double ma[3];
  double *major_axis = inclusion->orientation[id][0];
  double g = material->g;
  double v = material->v;
  double c = (1.0 + v) * g * e / (1.0 - v);
  double **offset = material->offset;
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
    TruncatedSphericalInclusionStress(rx, r, z1, z2, ma, c, inside,
                                      local_stress);
  }

  INCSTR_AddTensor(local_stress, stress);
}