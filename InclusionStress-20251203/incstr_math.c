#include <math.h>
#include <stdio.h>

void INCSTR_NormalizeVector(double *a) {
  double v = 0.0;

  for (int i = 0; i < 3; i++) {
    v += a[i] * a[i];
  }
  v = sqrt(v);
  for (int i = 0; i < 3; i++) {
    a[i] /= v;
  }
}

void INCSTR_NormalizeTensor(double **a) {
  for (int i = 0; i < 3; i++) {
    INCSTR_NormalizeVector(a[i]);
  }
}

double INCSTR_VectorLength(int n, double a[3]) {
  double v = 0.0;

  for (int i = 0; i < n; i++) {
    v += a[i] * a[i];
  }
  return sqrt(v);
}

void INCSTR_VectorProduct(double a[3], double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void INCSTR_RotateCoordinateVector(double rotation_tensor[3][3], double coord[3], double rotated_coord[3]) {
  for (int i = 0; i < 3; i++) {
    double v = 0.0;
    for (int j = 0; j < 3; j++) {
      v += rotation_tensor[i][j] * coord[j];
    }
    rotated_coord[i] = v;
  }
}

void INCSTR_RotateStressTensor(double rotation_tensor[3][3], double stress[3][3],
                             double rotated_stress[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double v = 0.0;
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          v += rotation_tensor[k][i] * rotation_tensor[l][j] *
               stress[k][l];
        }
      }
      rotated_stress[i][j] += v;
    }
  }
}

void INCSTR_AddTensor(double a[3][3], double b[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      b[i][j] += a[i][j];
    }
  }
}