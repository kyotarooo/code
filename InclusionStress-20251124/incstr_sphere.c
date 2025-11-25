//=======================================================================
// Reference:
//-----------------------------------------------------------------------
// Authors:
// A.L. Kolenikova, M.Yu. Gutkin, A.E. Romanov
// Title:
// Analytical elstic models of finite cylindrical and truncated spherical
// inclusions
// Journal: International Journal of Solids and Structures, Vol. 143,
// (2018), pp.59-72
//------------------------------------------------------------------------
// Authors:
// T. Mura
// Book:
// Micromechanics of Defects in Solids, (1987), Martinus Nijhoff, Boston
//========================================================================
#include <stdio.h>

#include "incstr_math.h"
#include "incstr_struct.h"

static void RotationTensor(double rx[3], double tensor[3][3]) {
  double vr[3] = {0.31412344, 0.786871414, 0.9234198124};  // Random direction
  double a[3], ar;

  for (int i = 0; i < 3; i++) {
    tensor[0][i] = rx[i];
  }

  a[0] = tensor[0][1] * vr[2] - tensor[0][2] * vr[1];
  a[1] = tensor[0][2] * vr[0] - tensor[0][0] * vr[2];
  a[2] = tensor[0][0] * vr[1] - tensor[0][1] * vr[0];
  ar = INCSTR_VectorLength(3, a);
  for (int i = 0; i < 3; i++) {
    tensor[1][i] = a[i] / ar;
  }

  tensor[2][0] = tensor[0][1] * tensor[1][2] - tensor[0][2] * tensor[1][1];
  tensor[2][1] = tensor[0][2] * tensor[1][0] - tensor[0][0] * tensor[1][2];
  tensor[2][2] = tensor[0][0] * tensor[1][1] - tensor[0][1] * tensor[1][0];
}

static void SphericalInclusionStress(double x[3], double e, double r, double g,
                                     double v, int *inside,
                                     double stress[3][3]) {
  double sr = -1.0 * (4.0 * (1.0 + v) * g * e) / (3.0 * (1.0 - v));
  double st = (2.0 * (1.0 + v) * g * e) / (3.0 * (1.0 - v));
  double distance = INCSTR_VectorLength(3, x);
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double rx[3];
  double rotation_tensor[3][3];

  if (r > distance) {
    st *= -2.0;
    (*inside) = 1;
  } else {
    double rd = r / distance;
    double rd3 = rd * rd * rd;

    sr *= rd3;
    st *= rd3;
  }

  local_stress[0][0] = sr;
  local_stress[1][1] = st;
  local_stress[2][2] = st;
  if (distance < 1.0e-20) {
    rx[0] = 1.0;
    rx[1] = 0.0;
    rx[2] = 0.0;
  } else {
    for (int iDim = 0; iDim < 3; iDim++) {
      rx[iDim] = x[iDim] / distance;
    }
  }
  RotationTensor(rx, rotation_tensor);
  INCSTR_RotateStressTensor(rotation_tensor, local_stress, stress);
}

void INCSTR_SphericalInclusionStress(double x[3], int id,
                                     INCLUSION_t *inclusion,
                                     MATERIAL_t *material, int *inside,
                                     double stress[3][3]) {
  double e = inclusion->eigen_strain[id];
  double r = inclusion->size[id][0];
  double g = material->g;
  double v = material->v;
  double *px = inclusion->x[id];
  double **offset = material->offset;
  double local_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  for (int i = 0; i < material->nImages; i++) {
    double rx[3];
    for (int iDim = 0; iDim < 3; iDim++) {
      rx[iDim] = x[iDim] - (px[iDim] + offset[i][iDim]);
    }
    SphericalInclusionStress(rx, e, r, g, v, inside, local_stress);
  }

  INCSTR_AddTensor(local_stress, stress);
}