#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fem_shape.h"
#include "fem_struct.h"
#include "fem_surface.h"
#include "fem_util.h"

static void ComputeDistributedSurfaceLoad(int elementID, int surfaceID,
                                          int direction, double v, FEM_t *fem) {
  FEMELEMENT_t *element = &fem->element[elementID];
  double element_x[8][3];
  double gauss_u[3] = {-0.774596669241483377035853079957, 0.0,
                       0.774596669241483377035853079957};
  double gauss_w[3] = {0.555555555555555555555555555556,
                       0.888888888888888888888888888889,
                       0.555555555555555555555555555556};
  int node_index[6][8] = {
      {3, 2, 1, 0, 6, 5, 4, 7},     {0, 1, 13, 12, 4, 9, 16, 8},
      {1, 2, 14, 13, 5, 10, 17, 9}, {2, 3, 15, 14, 6, 11, 18, 10},
      {3, 0, 12, 15, 7, 8, 19, 11}, {12, 13, 14, 15, 16, 17, 18, 19}};

  for (int iNode = 0; iNode < 8; iNode++) {
    int nodeID = element->nodeID[node_index[surfaceID][iNode]];
    FEMNODE_t *node = &fem->node[nodeID];
    for (int iDim = 0; iDim < 3; iDim++) {
      element_x[iNode][iDim] = node->x[iDim];
    }
  }

  for (int ix = 0; ix < 3; ix++) {
    double ux = gauss_u[ix];
    double wx = gauss_w[ix];
    for (int iy = 0; iy < 3; iy++) {
      double uy = gauss_u[iy];
      double wy = gauss_w[iy];
      double shape[8];
      double dshape[8][2];
      double dx_du[3][2];
      double jacobian;
      double u[2] = {ux, uy};

      FEM_ShapeFunction2D(u, shape);
      FEM_DShapeFunction2D(u, dshape);

      for (int iDim = 0; iDim < 3; iDim++) {
        for (int jDim = 0; jDim < 2; jDim++) {
          double v = 0.0;
          for (int iNode = 0; iNode < 8; iNode++) {
            v += dshape[iNode][jDim] * element_x[iNode][iDim];
          }
          dx_du[iDim][jDim] = v;
        }
      }

      {
        double ax = dx_du[1][0] * dx_du[2][1] - dx_du[2][0] * dx_du[1][1];
        double ay = dx_du[2][0] * dx_du[0][1] - dx_du[0][0] * dx_du[2][1];
        double az = dx_du[0][0] * dx_du[1][1] - dx_du[1][0] * dx_du[0][1];

        jacobian = sqrt(ax * ax + ay * ay + az * az);
      }

      for (int iNode = 0; iNode < 8; iNode++) {
        int nodeID = element->nodeID[node_index[surfaceID][iNode]];

        fem->node[nodeID].bc_f[direction] +=
            shape[iNode] * v * wx * wy * jacobian;
      }
    }
  }
}

void FEM_ReadInput(FEM_t *fem, char *directory) {
  int surface_elementID;
  char file_name[256];
  FILE *fp;

  sprintf(file_name, "%s/mesh.inp", directory);
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    fprintf(stdout, "[LOG] No FEM data is found\n");
    fflush(stdout);
    fem->nNodes = fem->nElements = 0;
    return;
  }
  fprintf(stdout, "[LOG] Reading FEM inputs\n");
  fflush(stdout);

  fscanf(fp, "%d%d%*d%*d", &fem->nElements, &fem->nNodes);

  fem->element = (FEMELEMENT_t *)malloc(fem->nElements * sizeof(FEMELEMENT_t));
  fem->node = (FEMNODE_t *)malloc(fem->nNodes * sizeof(FEMNODE_t));

  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];

    element->nodeID = (int *)malloc(20 * sizeof(int));
    element->surface = (int *)malloc(6 * sizeof(int));
    element->f = (double *)malloc(60 * sizeof(double));
    element->f_grid = (double *)malloc(60 * sizeof(double));

    for (int iNode = 0; iNode < 20; iNode++) {
      int nodeID;

      fscanf(fp, "%d", &nodeID);
      element->nodeID[iNode] = nodeID - 1;
    }
  }

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    node->bc = (int *)malloc(3 * sizeof(int));
    node->x = (double *)malloc(3 * sizeof(double));
    node->u = (double *)malloc(3 * sizeof(double));
    node->bc_u = (double *)malloc(3 * sizeof(double));
    node->f = (double *)malloc(3 * sizeof(double));
    node->bc_f = (double *)malloc(3 * sizeof(double));
    node->stress = (double *)malloc(6 * sizeof(double));

    for (int iDim = 0; iDim < 3; iDim++) {
      node->bc[iDim] = 0;
      node->u[iDim] = 0.0;
      node->bc_u[iDim] = 0.0;
      node->f[iDim] = 0.0;
      node->bc_f[iDim] = 0.0;
    }
    for (int iDim = 0; iDim < 6; iDim++) {
      node->stress[iDim] = 0.0;
    }
  }

  fem->xmin = (double *)malloc(3 * sizeof(double));
  fem->xmax = (double *)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    fem->xmin[iDim] = 1.0e+20;
    fem->xmax[iDim] = -1.0e+20;

    for (int iNode = 0; iNode < fem->nNodes; iNode++) {
      double x;

      fscanf(fp, "%lf", &x);
      fem->node[iNode].x[iDim] = x;
      fem->xmin[iDim] = fmin(fem->xmin[iDim], x);
      fem->xmax[iDim] = fmax(fem->xmax[iDim], x);
    }
  }

  while (1) {
    int elementID, surfaceID, direction;
    double v;

    fscanf(fp, "%d", &elementID);

    if (elementID == -1) break;

    fscanf(fp, "%d%d%lf", &surfaceID, &direction, &v);
    elementID--;
    surfaceID--;
    direction--;
    ComputeDistributedSurfaceLoad(elementID, surfaceID, direction, v, fem);
  }

  while (1) {
    int nodeID, direction;
    double v;

    fscanf(fp, "%d", &nodeID);

    if (nodeID == -1) break;

    fscanf(fp, "%d%lf", &direction, &v);
    nodeID--;
    direction--;
    fem->node[nodeID].bc_f[direction] += v;
  }

  while (1) {
    int nodeID, direction;
    double v;

    fscanf(fp, "%d", &nodeID);

    if (nodeID == -1) break;

    fscanf(fp, "%d%lf", &direction, &v);
    nodeID--;
    direction--;
    fem->node[nodeID].bc[direction] = 1;
    fem->node[nodeID].bc_u[direction] = v;
  }

  fclose(fp);

  FEM_FindSurface(fem);

  fem->nSurface_elements = 0;
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];

    if (element->nSurfaces != 0) fem->nSurface_elements += 1;
  }

  fem->surface_element = (int *)malloc(fem->nSurface_elements * sizeof(int));
  fem->surface_surface = (int **)malloc(fem->nSurface_elements * sizeof(int *));
  surface_elementID = 0;
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    int surfaceID = 0;
    FEMELEMENT_t *element = &fem->element[iElement];

    if (element->nSurfaces == 0) continue;

    fem->surface_element[surface_elementID] = iElement;
    fem->surface_surface[surface_elementID] =
        (int *)malloc(element->nSurfaces * sizeof(int));

    for (int iSurface = 0; iSurface < 6; iSurface++) {
      if (element->surface[iSurface]) {
        fem->surface_surface[surface_elementID][surfaceID] = iSurface;
        surfaceID += 1;
      }
    }
    surface_elementID += 1;
  }
}
