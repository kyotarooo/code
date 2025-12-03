#include <stdio.h>

void DD_ShapeFunction1D(double u, double shape[2]) {
  // linear line element
  shape[0] = 0.5 * (1.0 - u);
  shape[1] = 0.5 * (1.0 + u);
}

void DD_ShapeFunction3D(double u[3], double shape[8]) {
  double ux = u[0];
  double uy = u[1];
  double uz = u[2];

  // linear hexahedral element
  shape[0] = 0.125 * (1.0 - ux) * (1.0 - uy) * (1.0 - uz);
  shape[1] = 0.125 * (1.0 + ux) * (1.0 - uy) * (1.0 - uz);
  shape[2] = 0.125 * (1.0 + ux) * (1.0 + uy) * (1.0 - uz);
  shape[3] = 0.125 * (1.0 - ux) * (1.0 + uy) * (1.0 - uz);
  shape[4] = 0.125 * (1.0 - ux) * (1.0 - uy) * (1.0 + uz);
  shape[5] = 0.125 * (1.0 + ux) * (1.0 - uy) * (1.0 + uz);
  shape[6] = 0.125 * (1.0 + ux) * (1.0 + uy) * (1.0 + uz);
  shape[7] = 0.125 * (1.0 - ux) * (1.0 + uy) * (1.0 + uz);
}
