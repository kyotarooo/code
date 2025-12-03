void FEM_ShapeFunction3D(double u[3], double shape[20]) {
  double ux = u[0];
  double uy = u[1];
  double uz = u[2];

  shape[0] =
      -0.125 * (1.0 - ux) * (1.0 - uy) * (1.0 - uz) * (2.0 + ux + uy + uz);
  shape[1] =
      -0.125 * (1.0 + ux) * (1.0 - uy) * (1.0 - uz) * (2.0 - ux + uy + uz);
  shape[2] =
      -0.125 * (1.0 + ux) * (1.0 + uy) * (1.0 - uz) * (2.0 - ux - uy + uz);
  shape[3] =
      -0.125 * (1.0 - ux) * (1.0 + uy) * (1.0 - uz) * (2.0 + ux - uy + uz);

  shape[4] = 0.25 * (1.0 - ux * ux) * (1.0 - uy) * (1.0 - uz);
  shape[5] = 0.25 * (1.0 + ux) * (1.0 - uy * uy) * (1.0 - uz);
  shape[6] = 0.25 * (1.0 - ux * ux) * (1.0 + uy) * (1.0 - uz);
  shape[7] = 0.25 * (1.0 - ux) * (1.0 - uy * uy) * (1.0 - uz);

  shape[8] = 0.25 * (1.0 - ux) * (1.0 - uy) * (1.0 - uz * uz);
  shape[9] = 0.25 * (1.0 + ux) * (1.0 - uy) * (1.0 - uz * uz);
  shape[10] = 0.25 * (1.0 + ux) * (1.0 + uy) * (1.0 - uz * uz);
  shape[11] = 0.25 * (1.0 - ux) * (1.0 + uy) * (1.0 - uz * uz);

  shape[12] =
      -0.125 * (1.0 - ux) * (1.0 - uy) * (1.0 + uz) * (2.0 + ux + uy - uz);
  shape[13] =
      -0.125 * (1.0 + ux) * (1.0 - uy) * (1.0 + uz) * (2.0 - ux + uy - uz);
  shape[14] =
      -0.125 * (1.0 + ux) * (1.0 + uy) * (1.0 + uz) * (2.0 - ux - uy - uz);
  shape[15] =
      -0.125 * (1.0 - ux) * (1.0 + uy) * (1.0 + uz) * (2.0 + ux - uy - uz);

  shape[16] = 0.25 * (1.0 - ux * ux) * (1.0 - uy) * (1.0 + uz);
  shape[17] = 0.25 * (1.0 + ux) * (1.0 - uy * uy) * (1.0 + uz);
  shape[18] = 0.25 * (1.0 - ux * ux) * (1.0 + uy) * (1.0 + uz);
  shape[19] = 0.25 * (1.0 - ux) * (1.0 - uy * uy) * (1.0 + uz);
}

void FEM_DShapeFunction3D(double u[3], double dshape[20][3]) {
  double ux = u[0];
  double uy = u[1];
  double uz = u[2];

  dshape[0][0] = 0.125 * (1.0 - uy) * (1.0 - uz) * (1.0 + 2.0 * ux + uy + uz);
  dshape[0][1] = 0.125 * (1.0 - uz) * (1.0 - ux) * (1.0 + ux + 2.0 * uy + uz);
  dshape[0][2] = 0.125 * (1.0 - ux) * (1.0 - uy) * (1.0 + ux + uy + 2.0 * uz);
  dshape[1][0] = 0.125 * (1.0 - uy) * (1.0 - uz) * (-1.0 + 2.0 * ux - uy - uz);
  dshape[1][1] = 0.125 * (1.0 - uz) * (1.0 + ux) * (1.0 - ux + 2.0 * uy + uz);
  dshape[1][2] = 0.125 * (1.0 + ux) * (1.0 - uy) * (1.0 - ux + uy + 2.0 * uz);
  dshape[2][0] = 0.125 * (1.0 + uy) * (1.0 - uz) * (-1.0 + 2.0 * ux + uy - uz);
  dshape[2][1] = 0.125 * (1.0 - uz) * (1.0 + ux) * (-1.0 + ux + 2.0 * uy - uz);
  dshape[2][2] = 0.125 * (1.0 + ux) * (1.0 + uy) * (1.0 - ux - uy + 2.0 * uz);
  dshape[3][0] = 0.125 * (1.0 + uy) * (1.0 - uz) * (1.0 + 2.0 * ux - uy + uz);
  dshape[3][1] = 0.125 * (1.0 - uz) * (1.0 - ux) * (-1.0 - ux + 2.0 * uy - uz);
  dshape[3][2] = 0.125 * (1.0 - ux) * (1.0 + uy) * (1.0 + ux - uy + 2.0 * uz);
  dshape[4][0] = -0.5 * ux * (1.0 - uy) * (1.0 - uz);
  dshape[4][1] = -0.25 * (1.0 - ux * ux) * (1.0 - uz);
  dshape[4][2] = -0.25 * (1.0 - ux * ux) * (1.0 - uy);
  dshape[5][0] = 0.25 * (1.0 - uy * uy) * (1.0 - uz);
  dshape[5][1] = -0.5 * uy * (1.0 + ux) * (1.0 - uz);
  dshape[5][2] = -0.25 * (1.0 + ux) * (1.0 - uy * uy);
  dshape[6][0] = -0.5 * ux * (1.0 + uy) * (1.0 - uz);
  dshape[6][1] = 0.25 * (1.0 - ux * ux) * (1.0 - uz);
  dshape[6][2] = -0.25 * (1.0 - ux * ux) * (1.0 + uy);
  dshape[7][0] = -0.25 * (1.0 - uy * uy) * (1.0 - uz);
  dshape[7][1] = -0.5 * uy * (1.0 - ux) * (1.0 - uz);
  dshape[7][2] = -0.25 * (1.0 - ux) * (1.0 - uy * uy);
  dshape[8][0] = -0.25 * (1.0 - uy) * (1.0 - uz * uz);
  dshape[8][1] = -0.25 * (1.0 - ux) * (1.0 - uz * uz);
  dshape[8][2] = -0.5 * uz * (1.0 - ux) * (1.0 - uy);
  dshape[9][0] = 0.25 * (1.0 - uy) * (1.0 - uz * uz);
  dshape[9][1] = -0.25 * (1.0 + ux) * (1.0 - uz * uz);
  dshape[9][2] = -0.5 * uz * (1.0 + ux) * (1.0 - uy);
  dshape[10][0] = 0.25 * (1.0 + uy) * (1.0 - uz * uz);
  dshape[10][1] = 0.25 * (1.0 + ux) * (1.0 - uz * uz);
  dshape[10][2] = -0.5 * uz * (1.0 + ux) * (1.0 + uy);
  dshape[11][0] = -0.25 * (1.0 + uy) * (1.0 - uz * uz);
  dshape[11][1] = 0.25 * (1.0 - ux) * (1.0 - uz * uz);
  dshape[11][2] = -0.5 * uz * (1.0 - ux) * (1.0 + uy);
  dshape[12][0] = 0.125 * (1.0 - uy) * (1.0 + uz) * (1.0 + 2.0 * ux + uy - uz);
  dshape[12][1] = 0.125 * (1.0 + uz) * (1.0 - ux) * (1.0 + ux + 2.0 * uy - uz);
  dshape[12][2] = 0.125 * (1.0 - ux) * (1.0 - uy) * (-1.0 - ux - uy + 2.0 * uz);
  dshape[13][0] = 0.125 * (1.0 - uy) * (1.0 + uz) * (-1.0 + 2.0 * ux - uy + uz);
  dshape[13][1] = 0.125 * (1.0 + uz) * (1.0 + ux) * (1.0 - ux + 2.0 * uy - uz);
  dshape[13][2] = 0.125 * (1.0 + ux) * (1.0 - uy) * (-1.0 + ux - uy + 2.0 * uz);
  dshape[14][0] = 0.125 * (1.0 + uy) * (1.0 + uz) * (-1.0 + 2.0 * ux + uy + uz);
  dshape[14][1] = 0.125 * (1.0 + uz) * (1.0 + ux) * (-1.0 + ux + 2.0 * uy + uz);
  dshape[14][2] = 0.125 * (1.0 + ux) * (1.0 + uy) * (-1.0 + ux + uy + 2.0 * uz);
  dshape[15][0] = 0.125 * (1.0 + uy) * (1.0 + uz) * (1.0 + 2.0 * ux - uy - uz);
  dshape[15][1] = 0.125 * (1.0 + uz) * (1.0 - ux) * (-1.0 - ux + 2.0 * uy + uz);
  dshape[15][2] = 0.125 * (1.0 - ux) * (1.0 + uy) * (-1.0 - ux + uy + 2.0 * uz);
  dshape[16][0] = -0.5 * ux * (1.0 - uy) * (1.0 + uz);
  dshape[16][1] = -0.25 * (1.0 - ux * ux) * (1.0 + uz);
  dshape[16][2] = 0.25 * (1.0 - ux * ux) * (1.0 - uy);
  dshape[17][0] = 0.25 * (1.0 - uy * uy) * (1.0 + uz);
  dshape[17][1] = -0.5 * uy * (1.0 + ux) * (1.0 + uz);
  dshape[17][2] = 0.25 * (1.0 + ux) * (1.0 - uy * uy);
  dshape[18][0] = -0.5 * ux * (1.0 + uy) * (1.0 + uz);
  dshape[18][1] = 0.25 * (1.0 - ux * ux) * (1.0 + uz);
  dshape[18][2] = 0.25 * (1.0 - ux * ux) * (1.0 + uy);
  dshape[19][0] = -0.25 * (1.0 - uy * uy) * (1.0 + uz);
  dshape[19][1] = -0.5 * uy * (1.0 - ux) * (1.0 + uz);
  dshape[19][2] = 0.25 * (1.0 - ux) * (1.0 - uy * uy);
}

void FEM_ShapeFunction2D(double u[2], double shape[8]) {
  double ux = u[0];
  double uy = u[1];
  double a1 = 1.0 + ux;
  double a2 = 1.0 - ux;
  double b1 = 1.0 + uy;
  double b2 = 1.0 - uy;
  double c1 = 1.0 + ux + uy;
  double c2 = 1.0 + ux - uy;
  double c3 = 1.0 - ux + uy;
  double c4 = 1.0 - ux - uy;
  double d1 = 1.0 - ux * ux;
  double d2 = 1.0 - uy * uy;

  shape[0] = -0.25 * a2 * b2 * c1;
  shape[1] = -0.25 * a1 * b2 * c3;
  shape[2] = -0.25 * a1 * b1 * c4;
  shape[3] = -0.25 * a2 * b1 * c2;
  shape[4] = 0.5 * d1 * b2;
  shape[5] = 0.5 * d2 * a1;
  shape[6] = 0.5 * d1 * b1;
  shape[7] = 0.5 * d2 * a2;
}

void FEM_DShapeFunction2D(double u[2], double dshape[8][2]) {
  double ux = u[0];
  double uy = u[1];
  double a1 = 1.0 + ux;
  double a2 = 1.0 - ux;
  double b1 = 1.0 + uy;
  double b2 = 1.0 - uy;
  double c1 = 1.0 + ux + uy;
  double c2 = 1.0 + ux - uy;
  double c3 = 1.0 - ux + uy;
  double c4 = 1.0 - ux - uy;
  double d1 = 1.0 - ux * ux;
  double d2 = 1.0 - uy * uy;

  dshape[0][0] = -0.25 * (a2 * b2 - b2 * c1);
  dshape[0][1] = -0.25 * (a2 * b2 - a2 * c1);
  dshape[1][0] = -0.25 * (b2 * c3 - a1 * b2);
  dshape[1][1] = -0.25 * (a1 * b2 - a1 * c3);
  dshape[2][0] = -0.25 * (b1 * c4 - a1 * b1);
  dshape[2][1] = -0.25 * (a1 * c4 - a1 * b1);
  dshape[3][0] = -0.25 * (a2 * b1 - b1 * c2);
  dshape[3][1] = -0.25 * (a2 * c2 - a2 * b1);
  dshape[4][0] = -1.0 * ux * b2;
  dshape[4][1] = -0.5 * d1;
  dshape[5][0] = 0.5 * d2;
  dshape[5][1] = -1.0 * uy * a1;
  dshape[6][0] = -1.0 * ux * b1;
  dshape[6][1] = 0.5 * d1;
  dshape[7][0] = -0.5 * d2;
  dshape[7][1] = -1.0 * uy * a2;
}
