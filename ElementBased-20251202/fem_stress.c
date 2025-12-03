#include <math.h>
#include <stdio.h>

#include "fem_b_matrix.h"
#include "fem_d_matrix.h"
#include "fem_shape.h"
#include "fem_struct.h"
#include "fem_util.h"

void FEM_ComputeStress(FEM_t *fem) {
  double gauss_u[3] = {-0.7745966692414833770359, 0.0,
                       0.7745966692414833770359};
  double d_matrix[6][6];
  int node_integral_index[27] = {0, 8,  12, 7,  -1, 19, 3, 11, 15,
                                 4, -1, 16, -1, -1, -1, 6, -1, 18,
                                 1, 9,  13, 5,  -1, 17, 2, 10, 14};

  if (fem->nElements == 0) return;

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    node->iTmp = 0;
    for (int iDim = 0; iDim < 6; iDim++) {
      node->stress[iDim] = 0.0;
    }
  }

  FEM_MakeDMatrix(fem->material.e, fem->material.v, d_matrix);

  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];
    double element_x[20][3], element_u[20][3];
    int integralID = 0;

    for (int iNode = 0; iNode < 20; iNode++) {
      int nodeID = element->nodeID[iNode];
      FEMNODE_t *node = &fem->node[nodeID];
      double *x = node->x;
      double *u = node->u;
      for (int iDim = 0; iDim < 3; iDim++) {
        element_x[iNode][iDim] = x[iDim];
        element_u[iNode][iDim] = u[iDim];
      }
    }

    for (int ix = 0; ix < 3; ix++) {
      double ux = gauss_u[ix];

      for (int iy = 0; iy < 3; iy++) {
        double uy = gauss_u[iy];

        for (int iz = 0; iz < 3; iz++) {
          double uz = gauss_u[iz];
          double b_matrix[6][60];
          double jacobian;
          double strain[6];
          int node_index = node_integral_index[integralID];
          int integral_nodeID;
          FEMNODE_t *integral_node;
          double u[3] = {ux, uy, uz};

          integralID++;
          if (node_index == -1) continue;
          integral_nodeID = element->nodeID[node_index];
          integral_node = &fem->node[integral_nodeID];

          FEM_MakeBMatrix(b_matrix, &jacobian, u, element_x);

          for (int i = 0; i < 6; i++) {
            double v = 0.0;

            for (int iNode = 0; iNode < 20; iNode++) {
              for (int iDim = 0; iDim < 3; iDim++) {
                v += b_matrix[i][3 * iNode + iDim] * element_u[iNode][iDim];
              }
            }
            strain[i] = v;
          }

          for (int i = 0; i < 6; i++) {
            double v = 0.0;

            for (int j = 0; j < 6; j++) {
              v += d_matrix[i][j] * strain[j];
            }
            integral_node->stress[i] += v;
          }
          integral_node->iTmp++;
        }
      }
    }
  }

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    for (int iDim = 0; iDim < 6; iDim++) {
      node->stress[iDim] /= (double)node->iTmp;
    }
  }
}

static int FindElement(double x[3], double gx[3], FEM_t *fem) {
  int nearest_elementID = -1;
  double rmin = 1.0e+20;

  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    int found = 1;
    FEMELEMENT_t *element = &fem->element[iElement];

    for (int iSurface = 0; iSurface < 6; iSurface++) {
      double inner_product = 0.0;

      for (int iDim = 0; iDim < 3; iDim++) {
        double rx = element->surface_x[iSurface][iDim] - x[iDim];

        inner_product += rx * element->surface_nx[iSurface][iDim];
      }

      if (inner_product < 0.0) {
        found = 0;
        break;
      }
    }

    if (found) {
      FEM_NewtonRaphsonMethod(x, gx, element, fem);
      return iElement;
    }

    {
      double r = 0.0;

      for (int iDim = 0; iDim < 3; iDim++) {
        double rx = element->x[iDim] - x[iDim];

        r += rx * rx;
      }

      if (rmin > r) {
        rmin = r;
        nearest_elementID = iElement;
      }
    }
  }

  // どの要素にも含まれない場合は最も近い要素の番号を返す
  // 積分点はとりあえず原点とした
  gx[0] = 0.0;
  gx[1] = 0.0;
  gx[2] = 0.0;

  return nearest_elementID;
}

void FEM_CorrectionStress(double x[3], FEM_t *fem, double stress[3][3]) {
  double gx[3];
  int elementID;

  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      stress[iDim][jDim] = 0.0;
    }
  }

  if (fem->nElements == 0) return;

  elementID = FindElement(x, gx, fem);
  if (elementID != -1) {
    double shape[20];
    FEMELEMENT_t *element = &fem->element[elementID];
    double v00 = 0.0;
    double v11 = 0.0;
    double v22 = 0.0;
    double v01 = 0.0;
    double v12 = 0.0;
    double v02 = 0.0;

    FEM_ShapeFunction3D(gx, shape);
    for (int iNode = 0; iNode < 20; iNode++) {
      int nodeID = element->nodeID[iNode];
      FEMNODE_t *node = &fem->node[nodeID];
      double *node_stress = node->stress;
      double s = shape[iNode];

      v00 += s * node_stress[0];
      v11 += s * node_stress[1];
      v22 += s * node_stress[2];
      v01 += s * node_stress[3];
      v12 += s * node_stress[4];
      v02 += s * node_stress[5];
    }
    stress[0][0] = v00;
    stress[1][1] = v11;
    stress[2][2] = v22;
    stress[0][1] = stress[1][0] = v01;
    stress[1][2] = stress[2][1] = v12;
    stress[0][2] = stress[2][0] = v02;
  }
}
