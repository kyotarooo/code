#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_defs.h"
#include "dd_periodic_bc.h"
#include "dd_shape.h"
#include "dd_struct.h"

static void InitializeTensor(double a[3][3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      a[iDim][jDim] = 0.0;
    }
  }
}

static void SubstituteVector(double a[3], double b[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] = a[iDim];
  }
}

static void InterpolateStress(int gridID[3], double u[3], int *nx,
                              GRID_t ***grid, double stress[3][3]) {
  double shape[8];
  int grid_index[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                          {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

  // Shape function of linear hexahedral element
  DD_ShapeFunction3D(u, shape);

  // Interplate stress at 8 grid points
  for (int iGrid = 0; iGrid < 8; iGrid++) {
    double s = shape[iGrid];
    int target_gridID[3];
    GRID_t *target_grid;

    for (int iDim = 0; iDim < 3; iDim++) {
      target_gridID[iDim] = gridID[iDim] + grid_index[iGrid][iDim];
    }
    target_grid = &grid[target_gridID[0]][target_gridID[1]][target_gridID[2]];
    for (int iDim = 0; iDim < 3; iDim++) {
      for (int jDim = 0; jDim < 3; jDim++) {
        stress[iDim][jDim] += s * target_grid->stress[iDim][jDim];
      }
    }
  }
}

static void Offset(int *nx, double *size, double x[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    int id = (int)floor(x[iDim] / size[iDim]);
    if (id < 0)
      id = 0;
    else if (id >= nx[iDim])
      id = nx[iDim] - 1;

    x[iDim] -= (double)id * size[iDim];
  }
}

static void GridCoordinate(double x[3], int *nx, double *size, int gridID[3],
                           double u[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    // As for grid ID
    gridID[iDim] = (int)floor(x[iDim] / size[iDim]);
    if (gridID[iDim] < 0)
      gridID[iDim] = 0;
    else if (gridID[iDim] >= nx[iDim] - 1)
      gridID[iDim] = nx[iDim] - 2;

    // Regularized coordinate in element
    u[iDim] = 2.0 * (x[iDim] / size[iDim] - (double)gridID[iDim]) - 1.0;
    if (u[iDim] < -1.0)
      u[iDim] = -1.0;
    else if (u[iDim] > 1.0)
      u[iDim] = 1.0;
  }
}

static void SubgridStress(int depth, double ui[3], GRID_CELL_t *grid_cell, GRID_t ***grid,
                          double stress[3][3]) {
  int *nx = grid_cell->nx[depth];
  double size[3];
  double x[3], u[3];
  int gridID[3];
  GRID_t *igrid;

  // ui is ranging from -1 to 1.
  // x must be ranging from 0 to 2.
  for (int iDim = 0; iDim < 3; iDim++) {
    x[iDim] = ui[iDim] + 1.0;
  }

  // Grid size
  for(int iDim = 0; iDim < 3; iDim++) {
    size[iDim] = 2.0 / (double)(nx[iDim] - 1);
  }

  // Grid ID and regularized coordinate
  GridCoordinate(x, nx, size, gridID, u);
  igrid = &grid[gridID[0]][gridID[1]][gridID[2]];

  // Interpolate the stress using the shape function
  if (igrid->subgrid == NULL) {
    InterpolateStress(gridID, u, nx, grid, stress);
  } else {
    SubgridStress(depth + 1, u, grid_cell, igrid->subgrid, stress);
  }
}

void DD_GridStress(double xi[3], DD_t *dd, double stress[3][3]) {
  int depth = 0;
  int gridID[3];
  double x[3], u[3];
  INITIAL_STRESS_t *initial_stress = &dd->initial_stress;
  GRID_CELL_t *grid_cell = &initial_stress->grid_cell;
  int *nx;
  double size[3];
  GRID_t ***grid = grid_cell->grid;
  GRID_t *igrid;

  InitializeTensor(stress);

  // If no grid, nothing is done.
  if (!initial_stress->have_grids) return;

  nx = grid_cell->nx[depth];
  for(int iDim = 0; iDim < 3; iDim++) {
    size[iDim] = grid_cell->size[iDim];
  }

  // Integral point must be within the simulation volume
  SubstituteVector(xi, x);
  DD_PeriodicBoundaryConditionPosition(x, dd);

  // Grid ID and regularized coordinate
  Offset(initial_stress->nx, initial_stress->size, x);
  GridCoordinate(x, nx, size, gridID, u);
  igrid = &grid[gridID[0]][gridID[1]][gridID[2]];

  if (igrid->subgrid == NULL) {
    InterpolateStress(gridID, u, nx, grid, stress);
  } else {
    SubgridStress(depth + 1, u, grid_cell, igrid->subgrid, stress);
  }
}
