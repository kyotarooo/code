#include <stdio.h>

static void MatrixVectorProduct(double a[3][3], double b[3], double c[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    double v = 0.0;
    for (int jDim = 0; jDim < 3; jDim++) {
      v += a[iDim][jDim] * b[jDim];
    }
    c[iDim] = v;
  }
}

static void VectorProduct(double a[3], double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void DD_PeachKoehlerForce(double stress[3][3], double b[3], double t[3],
                          double f[3]) {
  double sb[3];

  MatrixVectorProduct(stress, b, sb);
  VectorProduct(sb, t, f);
}
