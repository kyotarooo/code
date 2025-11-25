#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "incstr_axial_cylinder.h"
#include "incstr_cuboid.h"
#include "incstr_cylinder.h"
#include "incstr_defs.h"
#include "incstr_dilatation_center.h"
#include "incstr_elastic_dipole.h"
#include "incstr_grid.h"
#include "incstr_sphere.h"
#include "incstr_struct.h"
#include "incstr_truncated_sphere.h"

static void InclusionStress(double x[3], INCLUSION_t *inclusion,
                            MATERIAL_t *material, int *inside,
                            double stress[3][3]) {
  for (int i = 0; i < inclusion->n; i++) {
    int shape = inclusion->shape[i];

    switch (shape) {
      case SPHERE:
        INCSTR_SphericalInclusionStress(x, i, inclusion, material, inside,
                                        stress);
        break;
      case CYLINDER:
        INCSTR_CylindricalInclusionStress(x, i, inclusion, material, inside,
                                          stress);
        break;
      case AXIAL_CYLINDER:
        INCSTR_AxialCylindricalInclusionStress(x, i, inclusion, material,
                                               inside, stress);
        break;
      case CUBOID:
        INCSTR_CuboidalInclusionStress(x, i, inclusion, material, inside,
                                       stress);
        break;
      case TRUNCATED_SPHERE:
        INCSTR_TruncatedSphericalInclusionStress(x, i, inclusion, material,
                                                 inside, stress);
        break;
      case ELASTIC_DIPOLE:
        INCSTR_ElasticDipoleStress(x, i, inclusion, material, stress);
        break;
      case DILATATION_CENTER:
        INCSTR_DilatationCenterStress(x, i, inclusion, material, stress);
        break;
    }
  }
}

static void MakeSubgrid(int depth, GRID_CELL_t *grid_cell, GRID_t *grid) {
  int *nx = grid_cell->nx[depth];

  grid->subgrid = (GRID_t ***)malloc(nx[0] * sizeof(GRID_t **));
  for (int ix = 0; ix < nx[0]; ix++) {
    grid->subgrid[ix] = (GRID_t **)malloc(nx[1] * sizeof(GRID_t *));
    for (int iy = 0; iy < nx[1]; iy++) {
      grid->subgrid[ix][iy] = (GRID_t *)malloc(nx[2] * sizeof(GRID_t));
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *subgrid = &grid->subgrid[ix][iy][iz];

        subgrid->x = (double *)malloc(3 * sizeof(double));
        subgrid->stress = (double *)malloc(6 * sizeof(double));
        for (int iDim = 1; iDim < 3; iDim++) {
          subgrid->stress[iDim] = subgrid->stress[iDim - 1] + 3;
        }
        subgrid->subgrid = NULL;
      }
    }
  }
}

static void MakeSubgrids(int depth, double *x0, GRID_CELL_t *grid_cell,
                         INCLUSION_t *inclusion, GRID_t ***grid) {
  int *nx = grid_cell->nx[depth];
  double *size = grid_cell->size[depth];

  for (int ix = 0; ix < nx[0] - 1; ix++) {
    for (int iy = 0; iy < nx[1] - 1; iy++) {
      for (int iz = 0; iz < nx[2] - 1; iz++) {
        int sum_inside = 0;

        for (int jx = 0; jx < 2; jx++) {
          for (int jy = 0; jy < 2; jy++) {
            for (int jz = 0; jz < 2; jz++) {
              sum_inside += grid[ix + jx][iy + jy][iz + jz].inside;
            }
          }
        }

        if (sum_inside != 0 && sum_inside != 8) {
          MakeSubgrid(depth + 1, grid_cell, &grid[ix][iy][iz]);
        }
      }
    }
  }

// 点欠陥がどのグリッドにあるか判定　→ そこのグリッドのサブグリッドを作成
  for (int i = 0; i < inclusion->n; i++) {
    int shape = inclusion->shape[i];

    if (shape == ELASTIC_DIPOLE || shape == DILATATION_CENTER) {
      double *x = inclusion->x[i];
      int id[3];
      GRID_t *igrid;

      for (int iDim = 0; iDim < 3; iDim++) {
        id[iDim] = (int)floor((x[iDim] - x0[iDim]) / size[iDim]);

        if (id[iDim] < 0) {
          id[iDim] = 0;
        } else if (id[iDim] >= nx[iDim] - 1) {
          id[iDim] = nx[iDim] - 2;
        }
      }

      igrid = &grid[id[0]][id[1]][id[2]];
      if (igrid->subgrid == NULL) {
        MakeSubgrid(depth + 1, grid_cell, igrid);
      }
    }
  }
}

static void SubgridStress(int depth, double *x0, GRID_CELL_t *grid_cell,
                          INCLUSION_t *inclusion, MATERIAL_t *material,
                          GRID_t ***grid) {
  int *nx = grid_cell->nx[depth];
  double *size = grid_cell->size[depth];
  int index[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {2, 0}};
  double x[3];

  for (int ix = 0; ix < nx[0]; ix++) {
    x[0] = x0[0] + (double)ix * size[0];

    for (int iy = 0; iy < nx[1]; iy++) {
      x[1] = x0[1] + (double)iy * size[1];

      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid[ix][iy][iz];
        double stress[3][3] = {
            {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        int inside = 0;

        x[2] = x0[2] + (double)iz * size[2];

        for (int iDim = 0; iDim < 3; iDim++) {
          igrid->x[iDim] = x[iDim];
        }

        InclusionStress(x, inclusion, material, &inside, stress);
        igrid->inside = inside;
        for (int i = 0; i < 6; i++) {
          igrid->stress[i] = stress[index[i][0]][index[i][1]];
        }
      }
    }
  }

  if (grid_cell->max_depth == depth + 1) return;

  MakeSubgrids(depth, x0, grid_cell, inclusion, grid);

  for (int ix = 0; ix < nx[0]; ix++) {
    for (int iy = 0; iy < nx[1]; iy++) {
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid[ix][iy][iz];

        if (igrid->subgrid != NULL) {
          SubgridStress(depth + 1, igrid->x, grid_cell, inclusion, material,
                        igrid->subgrid);
        }
      }
    }
  }
}

void INCSTR_InclusionStress(INCLUSION_t *inclusion, MATERIAL_t *material,
                            GRID_CELL_t *grid_cell) {
  int depth = 0;
  int *nx = grid_cell->nx[depth];
  double *size = grid_cell->size[depth];
  int index[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {2, 0}};
  double x0[3] = {0.0, 0.0, 0.0};

  grid_cell->grid = (GRID_t ***)malloc(nx[0] * sizeof(GRID_t **));
  for (int ix = 0; ix < nx[0]; ix++) {
    grid_cell->grid[ix] = (GRID_t **)malloc(nx[1] * sizeof(GRID_t *));
    for (int iy = 0; iy < nx[1]; iy++) {
      grid_cell->grid[ix][iy] = (GRID_t *)malloc(nx[2] * sizeof(GRID_t));
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid_cell->grid[ix][iy][iz];

        igrid->x = (double *)malloc(3 * sizeof(double));
        igrid->stress = (double *)malloc(6 * sizeof(double));
        igrid->subgrid = NULL;
      }
    }
  }

#pragma omp parallel for
  for (int ix = 0; ix < nx[0]; ix++) {
    double xi = (double)ix * size[0];
    for (int iy = 0; iy < nx[1]; iy++) {
      double yi = (double)iy * size[1];
      for (int iz = 0; iz < nx[2]; iz++) {
        double zi = (double)iz * size[2];
        double x[3] = {xi, yi, zi};
        double stress[3][3] = {
            {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        GRID_t *igrid = &grid_cell->grid[ix][iy][iz];
        int inside = 0;

        for (int iDim = 0; iDim < 3; iDim++) {
          igrid->x[iDim] = x[iDim];
        }

        InclusionStress(x, inclusion, material, &inside, stress);
        igrid->inside = inside;
        for (int i = 0; i < 6; i++) {
          igrid->stress[i] = stress[index[i][0]][index[i][1]];
        }
      }
    }
  }

  if (grid_cell->max_depth == depth + 1) return;

  MakeSubgrids(depth, x0, grid_cell, inclusion, grid_cell->grid);

#pragma omp parallel for
  for (int ix = 0; ix < nx[0]; ix++) {
    for (int iy = 0; iy < nx[1]; iy++) {
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid_cell->grid[ix][iy][iz];

        if (igrid->subgrid != NULL) {
          SubgridStress(depth + 1, igrid->x, grid_cell, inclusion, material,
                        igrid->subgrid);
        }
      }
    }
  }
}