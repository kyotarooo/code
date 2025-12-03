#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_dislocation_stress.h"
#include "dd_grid_stress.h"
#include "dd_node_handle.h"
#include "dd_struct.h"
#include "fem_shape.h"
#include "fem_solve.h"
#include "fem_stress.h"
#include "fem_struct.h"

static double VectorLength(double a[3]) {
  double r = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    r += a[iDim] * a[iDim];
  }
  return sqrt(r);
}

static double Jacobian(double dx_du[3][2], double n[3]) {
  double jacobian = 0.0;

  n[0] = dx_du[1][0] * dx_du[2][1] - dx_du[2][0] * dx_du[1][1];
  n[1] = dx_du[2][0] * dx_du[0][1] - dx_du[0][0] * dx_du[2][1];
  n[2] = dx_du[0][0] * dx_du[1][1] - dx_du[1][0] * dx_du[0][1];

  jacobian = VectorLength(n);
  for (int iDim = 0; iDim < 3; iDim++) {
    n[iDim] /= jacobian;
  }
  return jacobian;
}

static void IntegralPointPosition(double shape[8], double element_x[8][3],
                                  double x[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    double v = 0.0;

    for (int iNode = 0; iNode < 8; iNode++) {
      v += shape[iNode] * element_x[iNode][iDim];
    }
    x[iDim] = v;
  }
}

static void DIntegralPointPosition(double dshape[8][2], double element_x[8][3],
                                   double dx_du[3][2]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 2; jDim++) {
      double v = 0.0;
      for (int iNode = 0; iNode < 8; iNode++) {
        v += dshape[iNode][jDim] * element_x[iNode][iDim];
      }
      dx_du[iDim][jDim] = v;
    }
  }
}

static void Traction(double stress[3][3], double n[3], double t[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    double v = 0.0;
    for (int jDim = 0; jDim < 3; jDim++) {
      v += stress[iDim][jDim] * n[jDim];
    }
    t[iDim] = v;
  }
}

void DD_FEM_SurfaceTractionGridStress(DD_t *dd, FEM_t *fem) {
  double gauss_u[3] = {-0.774596669241483377035853079957,
                       0.000000000000000000000000000000,
                       0.774596669241483377035853079957};
  double gauss_w[3] = {0.555555555555555555555555555556,
                       0.888888888888888888888888888889,
                       0.555555555555555555555555555556};
  int node_index[6][8] = {
      {3, 2, 1, 0, 6, 5, 4, 7},     {0, 1, 13, 12, 4, 9, 16, 8},
      {1, 2, 14, 13, 5, 10, 17, 9}, {2, 3, 15, 14, 6, 11, 18, 10},
      {3, 0, 12, 15, 7, 8, 19, 11}, {12, 13, 14, 15, 16, 17, 18, 19}};

  if (fem->nElements == 0) return;

  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];
    for (int i = 0; i < 60; i++) element->f_grid[i] = 0.0;
  }

  if (!dd->initial_stress.have_grids) return;

  fprintf(stdout, "[LOG] Computing surface traction due to grid stress.\n");
  fflush(stdout);

#pragma omp parallel for
  for (int iElement = 0; iElement < fem->nSurface_elements; iElement++) {
    int elementID = fem->surface_element[iElement];
    FEMELEMENT_t *element = &fem->element[elementID];
    double *ef = element->f_grid;
    for (int iSurface = 0; iSurface < element->nSurfaces; iSurface++) {
      int surfaceID = fem->surface_surface[iElement][iSurface];
      double element_x[8][3];
      for (int iNode = 0; iNode < 8; iNode++) {
        int nodeID = element->nodeID[node_index[surfaceID][iNode]];
        FEMNODE_t *node = &fem->node[nodeID];
        double *x = node->x;
        for (int iDim = 0; iDim < 3; iDim++) {
          element_x[iNode][iDim] = x[iDim];
        }
      }

      for (int ix = 0; ix < 3; ix++) {
        double ux = gauss_u[ix];
        double wx = gauss_w[ix];
        for (int iy = 0; iy < 3; iy++) {
          double uy = gauss_u[iy];
          double wy = gauss_w[iy];
          double u[2] = {ux, uy};
          double shape[8];
          double dshape[8][2];
          double x[3];
          double dx_du[3][2];
          double stress[3][3] = {
              {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
          double jacobian, n[3];
          double t[3];

          FEM_ShapeFunction2D(u, shape);
          FEM_DShapeFunction2D(u, dshape);
          IntegralPointPosition(shape, element_x, x);
          DIntegralPointPosition(dshape, element_x, dx_du);
          jacobian = Jacobian(dx_du, n);

          DD_GridStress(x, dd, stress);
          Traction(stress, n, t);

          jacobian *= wx * wy;
          for (int iNode = 0; iNode < 8; iNode++) {
            int index = node_index[surfaceID][iNode];
            double sjww = shape[iNode] * jacobian;
            for (int iDim = 0; iDim < 3; iDim++) {
              ef[3 * index + iDim] -= sjww * t[iDim];
            }
          }
        }
      }
    }
  }
}

void DD_FEM_CorrectionField(DD_t *dd, FEM_t *fem) {
  double gauss_u[3] = {-0.774596669241483377035853079957,
                       0.000000000000000000000000000000,
                       0.774596669241483377035853079957};
  double gauss_w[3] = {0.555555555555555555555555555556,
                       0.888888888888888888888888888889,
                       0.555555555555555555555555555556};
  int node_index[6][8] = {
      {3, 2, 1, 0, 6, 5, 4, 7},     {0, 1, 13, 12, 4, 9, 16, 8},
      {1, 2, 14, 13, 5, 10, 17, 9}, {2, 3, 15, 14, 6, 11, 18, 10},
      {3, 0, 12, 15, 7, 8, 19, 11}, {12, 13, 14, 15, 16, 17, 18, 19}};

  if (fem->nElements == 0) return;

  fprintf(stdout, "[LOG] Updating the correction field\n");
  fflush(stdout);

  // To find a termination
  DD_FindElementsAroundNode(dd);

  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];
    for (int i = 0; i < 60; i++) element->f[i] = element->f_grid[i];
  }

#pragma omp parallel for
  for (int iElement = 0; iElement < fem->nSurface_elements; iElement++) {
    int elementID = fem->surface_element[iElement];
    FEMELEMENT_t *element = &fem->element[elementID];
    double *ef = element->f;
    for (int iSurface = 0; iSurface < element->nSurfaces; iSurface++) {
      int surfaceID = fem->surface_surface[iElement][iSurface];
      double element_x[8][3];
      for (int iNode = 0; iNode < 8; iNode++) {
        int nodeID = element->nodeID[node_index[surfaceID][iNode]];
        FEMNODE_t *node = &fem->node[nodeID];
        double *x = node->x;
        for (int iDim = 0; iDim < 3; iDim++) {
          element_x[iNode][iDim] = x[iDim];
        }
      }

      for (int ix = 0; ix < 3; ix++) {
        double ux = gauss_u[ix];
        double wx = gauss_w[ix];
        for (int iy = 0; iy < 3; iy++) {
          double uy = gauss_u[iy];
          double wy = gauss_w[iy];
          double u[2] = {ux, uy};
          double shape[8];
          double dshape[8][2];
          double x[3];
          double dx_du[3][2];
          double stress[3][3] = {
              {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
          double jacobian, n[3];
          double t[3];

          FEM_ShapeFunction2D(u, shape);
          FEM_DShapeFunction2D(u, dshape);
          IntegralPointPosition(shape, element_x, x);
          DIntegralPointPosition(dshape, element_x, dx_du);
          jacobian = Jacobian(dx_du, n);

          DD_DislocationStress(x, dd, stress);
          Traction(stress, n, t);

          jacobian *= wx * wy;
          for (int iNode = 0; iNode < 8; iNode++) {
            int index = node_index[surfaceID][iNode];
            double sjww = shape[iNode] * jacobian;
            for (int iDim = 0; iDim < 3; iDim++) {
              ef[3 * index + iDim] -= sjww * t[iDim];
            }
          }
        }
      }
    }
  }

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];
    double *f = node->f;
    for (int iDim = 0; iDim < 3; iDim++) {
      f[iDim] = 0.0;
    }
  }
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];
    double *ef = element->f;
    for (int iNode = 0; iNode < 20; iNode++) {
      int nodeID = element->nodeID[iNode];
      FEMNODE_t *node = &fem->node[nodeID];
      double *f = node->f;
      for (int iDim = 0; iDim < 3; iDim++) {
        f[iDim] += ef[3 * iNode + iDim];
      }
    }
  }
  FEM_ComputeSolution(fem);
  FEM_ComputeStress(fem);
}
