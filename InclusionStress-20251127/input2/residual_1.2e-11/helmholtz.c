#include <math.h>
#include <stdio.h>

int main(int argc, char **argv) {
  int n = 1000;
  double rmax = 2.0e-09;
  double dr = rmax / ((double)n - 1.0);
  double v = 0.3;
  double g = 80.0e+09;
  double q = 1.0e-27;
  double a0 = 3.15e-10;
  double c1 = 0.4 * a0;
  double c2 = 0.2 * a0;
  double c = g * q / (2.0 * M_PI) * (1.0 + v) / (1.0 - v);
  double d[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  for (int i = 0; i < n; i++) {
    double r = (double)i * dr + 1.0e-12;
    double rx[3] = {r, 0.0, 0.0};
    double stress[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double e1 = exp(-1.0 * r / c1);
    double e2 = exp(-1.0 * r / c2);
    double f2 = 1.0 -
                1.0 / (c1 * c1 - c2 * c2) * (c1 * c1 * e1 - c2 * c2 * e2) -
                r / (c1 * c1 - c2 * c2) * (c1 * e1 - c2 * e2) -
                r * r / (3.0 * (c1 * c1 - c2 * c2)) * (e1 - e2);
    double gbh = 1.0 / (4.0 * M_PI * (c1 * c1 - c2 * c2) * r) * (e1 - e2);

    for (int iDim = 0; iDim < 3; iDim++) {
      for (int jDim = 0; jDim < 3; jDim++) {
        stress[iDim][jDim] =
            c * ((d[iDim][jDim] / (r * r * r) -
                  3.0 * rx[iDim] * rx[jDim] / (r * r * r * r * r)) *
                     f2 -
                 8.0 * M_PI / 3.0 * d[iDim][jDim] * gbh);
      }
    }
    printf ("%e %e %e %e %e %e %e %e\n", r, f2, stress[0][0], stress[1][1], stress[2][2], stress[0][1], stress[1][2], stress[2][0]);
  }
  return 0;
}