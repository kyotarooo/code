#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fem_b_matrix.h"
#include "fem_d_matrix.h"
#include "fem_struct.h"

static void MakeMatrixIndex(FEM_t *fem) {
  int nNodes = fem->nNodes;
  FEMEQUATION_t *equation = &fem->equation;

  equation->nNodes = (int *)malloc(nNodes * sizeof(int));
  equation->nodeID = (int **)malloc(nNodes * sizeof(int *));

  for (int iNode = 0; iNode < nNodes; iNode++) {
    equation->nNodes[iNode] = 1;
  }

  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];
    for (int iNode = 0; iNode < 20; iNode++) {
      int iNodeID = element->nodeID[iNode];
      for (int jNode = 0; jNode < 20; jNode++) {
        int jNodeID = element->nodeID[jNode];

        if (iNodeID >= jNodeID) continue;

        equation->nNodes[iNodeID]++;
      }
    }
  }

  for (int iNode = 0; iNode < nNodes; iNode++) {
    equation->nodeID[iNode] =
        (int *)malloc(equation->nNodes[iNode] * sizeof(int));
    equation->nodeID[iNode][0] = iNode;
    equation->nNodes[iNode] = 1;
  }

  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];
    for (int iNode = 0; iNode < 20; iNode++) {
      int iNodeID = element->nodeID[iNode];
      for (int jNode = 0; jNode < 20; jNode++) {
        int kNode;
        int jNodeID = element->nodeID[jNode];
        int found = 0;

        if (iNodeID >= jNodeID) continue;

        for (kNode = 1; kNode < equation->nNodes[iNodeID]; kNode++) {
          if (equation->nodeID[iNodeID][kNode] == jNodeID) {
            found = 1;
            break;
          } else if (equation->nodeID[iNodeID][kNode] > jNodeID) {
            break;
          }
        }

        if (found) continue;

        for (int lNode = equation->nNodes[iNodeID] - 1; lNode >= kNode;
             lNode--) {
          equation->nodeID[iNodeID][lNode + 1] =
              equation->nodeID[iNodeID][lNode];
        }
        equation->nodeID[iNodeID][kNode] = jNodeID;
        equation->nNodes[iNodeID]++;
      }
    }
  }

  fem->equation.n = 0;
  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    fem->equation.n += 9 * fem->equation.nNodes[iNode] - 3;
  }

#ifndef _NO_FEM_
#ifdef _PARDISO_
  int count = 0;

  equation->nDofs = (int *)malloc(3 * nNodes * sizeof(int));
  equation->dofID = (int **)malloc(3 * nNodes * sizeof(int *));
  equation->start_index = (int *)malloc((3 * nNodes + 1) * sizeof(int));

  for (int iNode = 0; iNode < nNodes; iNode++) {
    int nNodesi = equation->nNodes[iNode];

    equation->nDofs[3 * iNode + 0] = 3 * nNodesi;
    equation->nDofs[3 * iNode + 1] = 3 * nNodesi - 1;
    equation->nDofs[3 * iNode + 2] = 3 * nNodesi - 2;
  }

  equation->dofID[0] = (int *)malloc(equation->n * sizeof(int));
  for (int iDof = 1; iDof < 3 * nNodes; iDof++) {
    equation->dofID[iDof] =
        equation->dofID[iDof - 1] + equation->nDofs[iDof - 1];
  }
  for (int iNode = 0; iNode < nNodes; iNode++) {
    int nNodesi = equation->nNodes[iNode];

    equation->dofID[3 * iNode + 0][0] = 3 * iNode + 0;
    equation->dofID[3 * iNode + 0][1] = 3 * iNode + 1;
    equation->dofID[3 * iNode + 0][2] = 3 * iNode + 2;
    equation->dofID[3 * iNode + 1][0] = 3 * iNode + 1;
    equation->dofID[3 * iNode + 1][1] = 3 * iNode + 2;
    equation->dofID[3 * iNode + 2][0] = 3 * iNode + 2;

    for (int jNode = 1; jNode < nNodes; jNode++) {
      int nodeID = equation->nodeID[iNode][jNode];

      equation->dofID[3 * iNode + 0][3 * jNode + 0 - 0] = 3 * nodeID + 0;
      equation->dofID[3 * iNode + 0][3 * jNode + 1 - 0] = 3 * nodeID + 1;
      equation->dofID[3 * iNode + 0][3 * jNode + 2 - 0] = 3 * nodeID + 2;
      equation->dofID[3 * iNode + 1][3 * jNode + 0 - 1] = 3 * nodeID + 0;
      equation->dofID[3 * iNode + 1][3 * jNode + 1 - 1] = 3 * nodeID + 1;
      equation->dofID[3 * iNode + 1][3 * jNode + 2 - 1] = 3 * nodeID + 2;
      equation->dofID[3 * iNode + 2][3 * jNode + 0 - 2] = 3 * nodeID + 0;
      equation->dofID[3 * iNode + 2][3 * jNode + 1 - 2] = 3 * nodeID + 1;
      equation->dofID[3 * iNode + 2][3 * jNode + 2 - 2] = 3 * nodeID + 2;
    }
  }
  equation->start_index[0] = 0;
  for (int iDof = 1; iDof < 3 * nNodes; iDof++) {
    equation->start_index[iDof] =
        equation->start_index[iDof - 1] + equation->nDofs[iDof - 1];
  }

  equation->bc = (int *)malloc(3 * nNodes * sizeof(int));
  for (int iNode = 0; iNode < nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    for (int iDim = 0; iDim < 3; iDim++) {
      if (node->bc[iDim])
        equation->bc[3 * iNode + iDim] = 1;
      else
        equation->bc[3 * iNode + iDim] = 0;
    }
  }

  equation->a = (double *)malloc(equation->n * sizeof(double));
  equation->b = (double *)malloc(3 * nNodes * sizeof(double));
  equation->x = (double *)malloc(3 * nNodes * sizeof(double));
  equation->b0 = (double *)malloc(3 * nNodes * sizeof(double));

  equation->MKL_system_index = (MKL_INT *)malloc(equation->n * sizeof(MKL_INT));
  for (int iDof = 0; iDof < 3 * nNodes; iDof++) {
    for (int jDof = 0; jDof < 3 * nNodes; jDof++) {
      equation->MKL_system_index[count] =
          (MKL_INT)equation->dofID[iDof][jDof] + 1;
      count++;
    }
  }

  equation->MKL_start_index =
      (MKL_INT *)malloc((3 * nNodes + 1) * sizeof(MKL_INT));
  for (int iDof = 0; iDof < 3 * nNodes; iDof++) {
    equation->MKL_start_index[iDof] = (MKL_INT)equation->start_index[iDof] + 1;
  }
  equation->MKL_start_index[3 * nNodes] = (MKL_INT)equation->n + 1;
#else
  fem->equation.start_index = (int *)malloc(3 * fem->nNodes * sizeof(int));
  fem->equation.row = (int *)malloc(fem->equation.n * sizeof(int));
  fem->equation.col = (int *)malloc(fem->equation.n * sizeof(int));

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    int nNodesi = fem->equation.nNodes[iNode];

    fem->equation.row[3 * iNode + 0] = 3 * nNodesi - 0;
    fem->equation.row[3 * iNode + 1] = 3 * nNodesi - 1;
    fem->equation.row[3 * iNode + 2] = 3 * nNodesi - 2;
  }

  fem->equation.start_index[0] = 0;
  for (int iNode = 1; iNode < 3 * fem->nNodes; iNode++) {
    fem->equation.start_index[iNode] =
        fem->equation.start_index[iNode - 1] + fem->equation.row[iNode - 1];
  }

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    int nNodesi = fem->equation.nNodes[iNode];
    int start_index0 = fem->equation.start_index[3 * iNode + 0];
    int start_index1 = fem->equation.start_index[3 * iNode + 1];
    int start_index2 = fem->equation.start_index[3 * iNode + 2];

    fem->equation.row[start_index0 + 0] = 3 * iNode + 0;
    fem->equation.col[start_index0 + 0] = 3 * iNode + 0;
    fem->equation.row[start_index0 + 1] = 3 * iNode + 0;
    fem->equation.col[start_index0 + 1] = 3 * iNode + 1;
    fem->equation.row[start_index0 + 2] = 3 * iNode + 0;
    fem->equation.col[start_index0 + 2] = 3 * iNode + 2;
    fem->equation.row[start_index1 + 0] = 3 * iNode + 1;
    fem->equation.col[start_index1 + 0] = 3 * iNode + 1;
    fem->equation.row[start_index1 + 1] = 3 * iNode + 1;
    fem->equation.col[start_index1 + 1] = 3 * iNode + 2;
    fem->equation.row[start_index2 + 0] = 3 * iNode + 2;
    fem->equation.col[start_index2 + 0] = 3 * iNode + 2;

    for (int jNode = 1; jNode < nNodesi; jNode++) {
      int nodeID = fem->equation.nodeID[iNode][jNode];

      fem->equation.row[start_index0 + 3 * jNode + 0 - 0] = 3 * iNode + 0;
      fem->equation.col[start_index0 + 3 * jNode + 0 - 0] = 3 * nodeID + 0;
      fem->equation.row[start_index0 + 3 * jNode + 1 - 0] = 3 * iNode + 0;
      fem->equation.col[start_index0 + 3 * jNode + 1 - 0] = 3 * nodeID + 1;
      fem->equation.row[start_index0 + 3 * jNode + 2 - 0] = 3 * iNode + 0;
      fem->equation.col[start_index0 + 3 * jNode + 2 - 0] = 3 * nodeID + 2;
      fem->equation.row[start_index1 + 3 * jNode + 0 - 1] = 3 * iNode + 1;
      fem->equation.col[start_index1 + 3 * jNode + 0 - 1] = 3 * nodeID + 0;
      fem->equation.row[start_index1 + 3 * jNode + 1 - 1] = 3 * iNode + 1;
      fem->equation.col[start_index1 + 3 * jNode + 1 - 1] = 3 * nodeID + 1;
      fem->equation.row[start_index1 + 3 * jNode + 2 - 1] = 3 * iNode + 1;
      fem->equation.col[start_index1 + 3 * jNode + 2 - 1] = 3 * nodeID + 2;
      fem->equation.row[start_index2 + 3 * jNode + 0 - 2] = 3 * iNode + 2;
      fem->equation.col[start_index2 + 3 * jNode + 0 - 2] = 3 * nodeID + 0;
      fem->equation.row[start_index2 + 3 * jNode + 1 - 2] = 3 * iNode + 2;
      fem->equation.col[start_index2 + 3 * jNode + 1 - 2] = 3 * nodeID + 1;
      fem->equation.row[start_index2 + 3 * jNode + 2 - 2] = 3 * iNode + 2;
      fem->equation.col[start_index2 + 3 * jNode + 2 - 2] = 3 * nodeID + 2;
    }
  }

  fem->equation.bc = (double *)malloc(3 * fem->nNodes * sizeof(double));
  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    for (int iDim = 0; iDim < 3; iDim++) {
      if (fem->node[iNode].bc[iDim])
        fem->equation.bc[3 * iNode + iDim] = 0.0;
      else
        fem->equation.bc[3 * iNode + iDim] = 1.0;
    }
  }

  fem->equation.a = (double *)malloc(fem->equation.n * sizeof(double));
  fem->equation.a0 = (double *)malloc(fem->equation.n * sizeof(double));
  fem->equation.b = (double *)malloc(3 * fem->nNodes * sizeof(double));
  fem->equation.x = (double *)malloc(3 * fem->nNodes * sizeof(double));
#endif
#endif

  // Initial value = 1
  fem->equation.scale_u = 1.0;
  fem->equation.scale_f = 1.0;
}

static void ComputeElementStiffnessMatrix(FEMELEMENT_t *element,
                                          double element_k[60][60],
                                          FEM_t *fem) {
  double element_x[20][3];
  double gauss_u[3] = {-0.774596669241483377035853079957, 0.0,
                       0.774596669241483377035853079957};
  double gauss_w[3] = {0.555555555555555555555555555556,
                       0.888888888888888888888888888889,
                       0.555555555555555555555555555556};
  double d_matrix[6][6];

  for (int ix = 0; ix < 60; ix++) {
    for (int iy = 0; iy < 60; iy++) {
      element_k[ix][iy] = 0.0;
    }
  }

  for (int iNode = 0; iNode < 20; iNode++) {
    int nodeID = element->nodeID[iNode];
    FEMNODE_t *node = &fem->node[nodeID];
    for (int iDim = 0; iDim < 3; iDim++) {
      element_x[iNode][iDim] = node->x[iDim];
    }
  }

  FEM_MakeDMatrix(fem->material.e, fem->material.v, d_matrix);

  for (int ix = 0; ix < 3; ix++) {
    double ux = gauss_u[ix];
    double wx = gauss_w[ix];
    for (int iy = 0; iy < 3; iy++) {
      double uy = gauss_u[iy];
      double wy = gauss_w[iy];
      for (int iz = 0; iz < 3; iz++) {
        double uz = gauss_u[iz];
        double wz = gauss_w[iz];
        double b_matrix[6][60];
        double jacobian;
        double db_matrix[6][60];
        double u[3] = {ux, uy, uz};

        FEM_MakeBMatrix(b_matrix, &jacobian, u, element_x);

        for (int i = 0; i < 6; i++) {
          for (int j = 0; j < 60; j++) {
            double v = 0.0;

            for (int k = 0; k < 6; k++) {
              v += d_matrix[i][k] * b_matrix[k][j];
            }
            db_matrix[i][j] = v;
          }
        }

        for (int i = 0; i < 60; i++) {
          for (int j = 0; j < 60; j++) {
            double v = 0.0;

            for (int k = 0; k < 6; k++) {
              v += b_matrix[k][i] * db_matrix[k][j];
            }
            element_k[i][j] += v * wx * wy * wz * jacobian;
          }
        }
        fem->volume += wx * wy * wz * jacobian;
      }
    }
  }
}

#ifndef _NO_FEM_
#ifdef _PARDISO_
static void AssembleGlobalStiffnessMatrix(FEMELEMENT_t *element,
                                          double element_k[60][60],
                                          FEMEQUATION_t *equation) {
  for (int iNode = 0; iNode < 20; iNode++) {
    int iNodeID = element->nodeID[iNode];
    int iDofID0 = 3 * iNodeID + 0;
    int iDofID1 = 3 * iNodeID + 1;
    int iDofID2 = 3 * iNodeID + 2;
    int start_index0 = equation->start_index[iDofID0];
    int start_index1 = equation->start_index[iDofID1];
    int start_index2 = equation->start_index[iDofID2];
    for (int jNode = 0; jNode < 20; jNode++) {
      int jNodeID = element->nodeID[jNode];
      int jDofID0 = 3 * jNodeID + 0;

      if (iNodeID > jNodeID) continue;

      for (int iDof = 0; iDof < equation->nDofs[iDofID0]; iDof += 3) {
        if (equation->dofID[iDofID0][iDof] != jDofID0) continue;

        if (iNodeID == jNodeID) {
          equation->a[start_index0 + iDof + 0] +=
              element_k[3 * iNode + 0][3 * jNode + 0];
          equation->a[start_index0 + iDof + 1] +=
              element_k[3 * iNode + 0][3 * jNode + 1];
          equation->a[start_index0 + iDof + 2] +=
              element_k[3 * iNode + 0][3 * jNode + 2];
          equation->a[start_index1 + iDof + 0] +=
              element_k[3 * iNode + 1][3 * jNode + 1];
          equation->a[start_index1 + iDof + 1] +=
              element_k[3 * iNode + 1][3 * jNode + 2];
          equation->a[start_index2 + iDof + 0] +=
              element_k[3 * iNode + 2][3 * jNode + 2];
        } else {
          equation->a[start_index0 + iDof + 0 - 0] +=
              element_k[3 * iNode + 0][3 * jNode + 0];
          equation->a[start_index0 + iDof + 1 - 0] +=
              element_k[3 * iNode + 0][3 * jNode + 1];
          equation->a[start_index0 + iDof + 2 - 0] +=
              element_k[3 * iNode + 0][3 * jNode + 2];
          equation->a[start_index1 + iDof + 0 - 1] +=
              element_k[3 * iNode + 1][3 * jNode + 0];
          equation->a[start_index1 + iDof + 1 - 1] +=
              element_k[3 * iNode + 1][3 * jNode + 1];
          equation->a[start_index1 + iDof + 2 - 1] +=
              element_k[3 * iNode + 1][3 * jNode + 2];
          equation->a[start_index2 + iDof + 0 - 2] +=
              element_k[3 * iNode + 2][3 * jNode + 0];
          equation->a[start_index2 + iDof + 1 - 2] +=
              element_k[3 * iNode + 2][3 * jNode + 1];
          equation->a[start_index2 + iDof + 2 - 2] +=
              element_k[3 * iNode + 2][3 * jNode + 2];
        }
      }
    }
  }
}
#else
static void AssembleGlobalStiffnessMatrix(FEMELEMENT_t *element,
                                          double element_k[60][60],
                                          FEMEQUATION_t *equation) {
  for (int iNode = 0; iNode < 20; iNode++) {
    int iNodeID = element->nodeID[iNode];
    int iDofID0 = 3 * iNodeID + 0;
    int iDofID1 = 3 * iNodeID + 1;
    int iDofID2 = 3 * iNodeID + 2;
    int start_index0 = equation->start_index[iDofID0];
    int start_index1 = equation->start_index[iDofID1];
    int start_index2 = equation->start_index[iDofID2];
    for (int jNode = 0; jNode < 20; jNode++) {
      int jNodeID = element->nodeID[jNode];
      int jDofID0 = 3 * jNodeID;

      if (iNodeID > jNodeID) continue;

      for (int iDof = 0; iDof < equation->n; iDof++) {
        if (equation->row[start_index0 + iDof] == iDofID0 &&
            equation->col[start_index0 + iDof] == jDofID0) {
          if (iNodeID == jNodeID) {
            equation->a[start_index0 + iDof + 0] +=
                element_k[3 * iNode + 0][3 * jNode + 0];
            equation->a[start_index0 + iDof + 1] +=
                element_k[3 * iNode + 0][3 * jNode + 1];
            equation->a[start_index0 + iDof + 2] +=
                element_k[3 * iNode + 0][3 * jNode + 2];
            equation->a[start_index1 + iDof + 0] +=
                element_k[3 * iNode + 1][3 * jNode + 1];
            equation->a[start_index1 + iDof + 1] +=
                element_k[3 * iNode + 1][3 * jNode + 2];
            equation->a[start_index2 + iDof + 0] +=
                element_k[3 * iNode + 2][3 * jNode + 2];
          } else {
            equation->a[start_index0 + iDof + 0 - 0] +=
                element_k[3 * iNode + 0][3 * jNode + 0];
            equation->a[start_index0 + iDof + 1 - 0] +=
                element_k[3 * iNode + 0][3 * jNode + 1];
            equation->a[start_index0 + iDof + 2 - 0] +=
                element_k[3 * iNode + 0][3 * jNode + 2];
            equation->a[start_index1 + iDof + 0 - 1] +=
                element_k[3 * iNode + 1][3 * jNode + 0];
            equation->a[start_index1 + iDof + 1 - 1] +=
                element_k[3 * iNode + 1][3 * jNode + 1];
            equation->a[start_index1 + iDof + 2 - 1] +=
                element_k[3 * iNode + 1][3 * jNode + 2];
            equation->a[start_index2 + iDof + 0 - 2] +=
                element_k[3 * iNode + 2][3 * jNode + 0];
            equation->a[start_index2 + iDof + 1 - 2] +=
                element_k[3 * iNode + 2][3 * jNode + 1];
            equation->a[start_index2 + iDof + 2 - 2] +=
                element_k[3 * iNode + 2][3 * jNode + 2];
          }
          break;
        }
      }
    }
  }
}
#endif
#else
static void AssembleGlobalStiffnessMatrix(FEMELEMENT_t *element,
                                          double element_k[60][60],
                                          FEMEQUATION_t *equation) {}
#endif

void FEM_StiffnessMatrix(FEM_t *fem) {
  FEMEQUATION_t *equation = &fem->equation;

  if (fem->nElements == 0) return;

  fprintf(stdout, "[LOG] Computing the stiffness matrix\n");
  fflush(stdout);

  MakeMatrixIndex(fem);
  for (int i = 0; i < equation->n; i++) equation->a[i] = 0.0;

  fem->volume = 0.0;
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];
    double element_k[60][60];

    ComputeElementStiffnessMatrix(element, element_k, fem);
    AssembleGlobalStiffnessMatrix(element, element_k, equation);
  }

#ifndef _NO_FEM_
#ifndef _PARDISO_
  for (int i = 0; i < fem->equation.n; i++) {
    fem->equation.a0[i] = fem->equation.a[i];
  }
#endif
#endif
}
