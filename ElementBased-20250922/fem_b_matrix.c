#include "fem_shape.h"

void FEM_MakeBMatrix(double b_matrix[6][60], double *jacobian, double u[3],
                     double element_x[20][3]) {
  double dshape[20][3];
  double dx_du[3][3];
  double du_dx[3][3];
  double dN_dx[20][3];
  double det;

  // 形状関数の微分の計算
  FEM_DShapeFunction3D(u, dshape);

  // ここからヤコビアンと形状関数の空間微分を求める
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      double v = 0.0;

      for (int iNode = 0; iNode < 20; iNode++) {
        v += dshape[iNode][jDim] * element_x[iNode][iDim];
      }
      dx_du[iDim][jDim] = v;
    }
  }

  det = dx_du[0][0] * (dx_du[1][1] * dx_du[2][2] - dx_du[1][2] * dx_du[2][1]) +
        dx_du[0][1] * (dx_du[1][2] * dx_du[2][0] - dx_du[1][0] * dx_du[2][2]) +
        dx_du[0][2] * (dx_du[1][0] * dx_du[2][1] - dx_du[1][1] * dx_du[2][0]);
  (*jacobian) = det;
  du_dx[0][0] = (dx_du[1][1] * dx_du[2][2] - dx_du[1][2] * dx_du[2][1]) / det;
  du_dx[0][1] = (dx_du[0][2] * dx_du[2][1] - dx_du[0][1] * dx_du[2][2]) / det;
  du_dx[0][2] = (dx_du[0][1] * dx_du[1][2] - dx_du[0][2] * dx_du[1][1]) / det;
  du_dx[1][0] = (dx_du[1][2] * dx_du[2][0] - dx_du[1][0] * dx_du[2][2]) / det;
  du_dx[1][1] = (dx_du[0][0] * dx_du[2][2] - dx_du[0][2] * dx_du[2][0]) / det;
  du_dx[1][2] = (dx_du[0][2] * dx_du[1][0] - dx_du[0][0] * dx_du[1][2]) / det;
  du_dx[2][0] = (dx_du[1][0] * dx_du[2][1] - dx_du[1][1] * dx_du[2][0]) / det;
  du_dx[2][1] = (dx_du[0][1] * dx_du[2][0] - dx_du[0][0] * dx_du[2][1]) / det;
  du_dx[2][2] = (dx_du[0][0] * dx_du[1][1] - dx_du[0][1] * dx_du[1][0]) / det;

  for (int iNode = 0; iNode < 20; iNode++) {
    for (int iDim = 0; iDim < 3; iDim++) {
      double v = 0.0;

      for (int jDim = 0; jDim < 3; jDim++) {
        v += dshape[iNode][jDim] * du_dx[jDim][iDim];
      }
      dN_dx[iNode][iDim] = v;
    }
  }

  // やっと変位ーひずみ行列の作成
  for (int iNode = 0; iNode < 20; iNode++) {
    b_matrix[0][3 * iNode + 0] = dN_dx[iNode][0];
    b_matrix[0][3 * iNode + 1] = 0.0;
    b_matrix[0][3 * iNode + 2] = 0.0;
    b_matrix[1][3 * iNode + 0] = 0.0;
    b_matrix[1][3 * iNode + 1] = dN_dx[iNode][1];
    b_matrix[1][3 * iNode + 2] = 0.0;
    b_matrix[2][3 * iNode + 0] = 0.0;
    b_matrix[2][3 * iNode + 1] = 0.0;
    b_matrix[2][3 * iNode + 2] = dN_dx[iNode][2];
    b_matrix[3][3 * iNode + 0] = dN_dx[iNode][1];
    b_matrix[3][3 * iNode + 1] = dN_dx[iNode][0];
    b_matrix[3][3 * iNode + 2] = 0.0;
    b_matrix[4][3 * iNode + 0] = 0.0;
    b_matrix[4][3 * iNode + 1] = dN_dx[iNode][2];
    b_matrix[4][3 * iNode + 2] = dN_dx[iNode][1];
    b_matrix[5][3 * iNode + 0] = dN_dx[iNode][2];
    b_matrix[5][3 * iNode + 1] = 0.0;
    b_matrix[5][3 * iNode + 2] = dN_dx[iNode][0];
  }
}
