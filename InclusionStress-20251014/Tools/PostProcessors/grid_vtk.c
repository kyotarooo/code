#include <stdio.h>
#include <stdlib.h>

typedef struct __grid_t_ {
  int id;
  double *x;
  struct __grid_t_ ***subgrid;
} GRID_t;

typedef struct __grid_cell_t_ {
  int max_depth;
  int **nx;
  double **size;
  GRID_t ***grid;
} GRID_CELL_t;

static void ReadSubgrids(int *gridID, double x0[3], int depth,
                         GRID_CELL_t *grid_cell, GRID_t *grid, FILE *fp) {
  int *nx = grid_cell->nx[depth];
  double *size = grid_cell->size[depth];

  grid->subgrid = (GRID_t ***)malloc(nx[0] * sizeof(GRID_t **));
  for (int ix = 0; ix < nx[0]; ix++) {
    double xi = x0[0] + (double)ix * size[0];

    grid->subgrid[ix] = (GRID_t **)malloc(nx[1] * sizeof(GRID_t *));
    for (int iy = 0; iy < nx[1]; iy++) {
      double yi = x0[1] + (double)iy * size[1];

      grid->subgrid[ix][iy] = (GRID_t *)malloc(nx[2] * sizeof(GRID_t));
      for (int iz = 0; iz < nx[2]; iz++) {
        double zi = x0[2] + (double)iz * size[2];
        double x[3] = {xi, yi, zi};
        GRID_t *subgrid = &grid->subgrid[ix][iy][iz];
        double stress[9];
        int have_subgrid;

        subgrid->x = (double *)malloc(3 * sizeof(double));
        for (int iDim = 0; iDim < 3; iDim++) {
          subgrid->x[iDim] = x[iDim];
        }

        fread(stress, sizeof(double), 9, fp);
        fread(&have_subgrid, sizeof(int), 1, fp);

        subgrid->id = (*gridID);
        (*gridID) += 1;

        if (have_subgrid) {
          ReadSubgrids(gridID, x, depth + 1, grid_cell, subgrid, fp);
        } else {
          subgrid->subgrid = NULL;
        }
      }
    }
  }
}

static void ReadGrids(GRID_CELL_t *grid_cell, char *directory) {
  char file_name[256];
  FILE *fp;
  double size[3] = {1.0, 1.0, 1.0};
  int depth = 0;
  int *nx;
  double *grid_size;
  int gridID = 0;

  sprintf(file_name, "%s/grid_stress.inp", directory);
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot open %s file\n", file_name);
    exit(1);
  }

  fread(&grid_cell->max_depth, sizeof(int), 1, fp);
  grid_cell->nx = (int **)malloc(grid_cell->max_depth * sizeof(int *));
  for (int i = 0; i < grid_cell->max_depth; i++) {
    grid_cell->nx[i] = (int *)malloc(3 * sizeof(int));
    fread(grid_cell->nx[i], sizeof(int), 3, fp);
  }
  nx = grid_cell->nx[depth];

  grid_cell->size = (double **)malloc(grid_cell->max_depth * sizeof(double *));
  for (int i = 0; i < grid_cell->max_depth; i++) {
    grid_cell->size[i] = (double *)malloc(3 * sizeof(double));
    for (int iDim = 0; iDim < 3; iDim++) {
      grid_cell->size[i][iDim] = size[iDim] / (double)(grid_cell->nx[i][iDim] - 1);
      size[iDim] = grid_cell->size[i][iDim];
    }
  }
  grid_size = grid_cell->size[depth];

  grid_cell->grid = (GRID_t ***)malloc(nx[0] * sizeof(GRID_t **));
  for (int ix = 0; ix < nx[0]; ix++) {
    double xi = (double)ix * grid_size[0];

    grid_cell->grid[ix] = (GRID_t **)malloc(nx[1] * sizeof(GRID_t *));
    for (int iy = 0; iy < nx[1]; iy++) {
      double yi = (double)iy * grid_size[1];

      grid_cell->grid[ix][iy] = (GRID_t *)malloc(nx[2] * sizeof(GRID_t));
      for (int iz = 0; iz < nx[2]; iz++) {
        double zi = (double)iz * grid_size[2];
        double x[3] = {xi, yi, zi};
        GRID_t *igrid = &grid_cell->grid[ix][iy][iz];
        double stress[9];
        int have_subgrid;

        igrid->x = (double *)malloc(3 * sizeof(double));
        for (int iDim = 0; iDim < 3; iDim++) {
          igrid->x[iDim] = x[iDim];
        }

        fread(stress, sizeof(double), 9, fp);
        fread(&have_subgrid, sizeof(int), 1, fp);

        igrid->id = gridID;
        gridID += 1;

        if (have_subgrid) {
          ReadSubgrids(&gridID, x, depth + 1, grid_cell, igrid, fp);
        } else {
          igrid->subgrid = NULL;
        }
      }
    }
  }

  fclose(fp);
}

static void NumberOfPointsInSubgrids(int *n, int depth, GRID_CELL_t *grid_cell,
                                    GRID_t ***grid) {
  int *nx = grid_cell->nx[depth];

  (*n) += nx[0] * nx[1] * nx[2];
  for (int ix = 0; ix < nx[0]; ix++) {
    for (int iy = 0; iy < nx[1]; iy++) {
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid[ix][iy][iz];

        if (igrid->subgrid != NULL) {
          NumberOfPointsInSubgrids(n, depth + 1, grid_cell, igrid->subgrid);
        }
      }
    }
  }
}

static int NumberOfPoints(GRID_CELL_t *grid_cell) {
  int depth = 0;
  int *nx = grid_cell->nx[depth];
  int n = 0;

  n = nx[0] * nx[1] * nx[2];
  for (int ix = 0; ix < nx[0]; ix++) {
    for (int iy = 0; iy < nx[1]; iy++) {
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid_cell->grid[ix][iy][iz];

        if (igrid->subgrid != NULL) {
          NumberOfPointsInSubgrids(&n, depth + 1, grid_cell, igrid->subgrid);
        }
      }
    }
  }

  return n;
}

static void NumberOfCellsInSubgrids(int *n, int depth, GRID_CELL_t *grid_cell,
                                    GRID_t ***grid) {
  int *nx = grid_cell->nx[depth];

  (*n) += (nx[0] - 1) * (nx[1] - 1) * (nx[2] - 1);
  for (int ix = 0; ix < nx[0]; ix++) {
    for (int iy = 0; iy < nx[1]; iy++) {
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid[ix][iy][iz];

        if (igrid->subgrid != NULL) {
          NumberOfCellsInSubgrids(n, depth + 1, grid_cell, igrid->subgrid);
        }
      }
    }
  }
}

static int NumberOfCells(GRID_CELL_t *grid_cell) {
  int depth = 0;
  int *nx = grid_cell->nx[depth];
  int n = 0;

  n = (nx[0] - 1) * (nx[1] - 1) * (nx[2] - 1);
  for (int ix = 0; ix < nx[0]; ix++) {
    for (int iy = 0; iy < nx[1]; iy++) {
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid_cell->grid[ix][iy][iz];

        if (igrid->subgrid != NULL) {
          NumberOfCellsInSubgrids(&n, depth + 1, grid_cell, igrid->subgrid);
        }
      }
    }
  }

  return n;
}

static void WriteSubgridPoints(int depth, GRID_CELL_t *grid_cell,
                               GRID_t ***grid, FILE *fp) {
  int *nx = grid_cell->nx[depth];

  for (int ix = 0; ix < nx[0]; ix++) {
    for (int iy = 0; iy < nx[1]; iy++) {
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid[ix][iy][iz];

        for (int iDim = 0; iDim < 3; iDim++) {
          fprintf(fp, "%e ", igrid->x[iDim]);
        }
        fprintf(fp, "\n");

        if (igrid->subgrid != NULL) {
          WriteSubgridPoints(depth + 1, grid_cell, igrid->subgrid, fp);
        }
      }
    }
  }
}

static void WriteSubgridCells(int depth, GRID_CELL_t *grid_cell, GRID_t ***grid,
                              FILE *fp) {
  int *nx = grid_cell->nx[depth];
  int idx[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                   {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

  for (int ix = 0; ix < nx[0] - 1; ix++) {
    for (int iy = 0; iy < nx[1] - 1; iy++) {
      for (int iz = 0; iz < nx[2] - 1; iz++) {
        GRID_t *igrid = &grid[ix][iy][iz];

        fprintf(fp, "8 ");
        for (int i = 0; i < 8; i++) {
          fprintf(fp, "%d ",
                  grid[ix + idx[i][0]][iy + idx[i][1]][iz + idx[i][2]].id);
        }

        if (igrid->subgrid != NULL) {
          WriteSubgridCells(depth + 1, grid_cell, igrid->subgrid, fp);
        }
      }
    }
  }
}

static void WriteGrids(GRID_CELL_t *grid_cell, char *directory) {
  int nPoints = NumberOfPoints(grid_cell);
  int nCells = NumberOfCells(grid_cell);
  char file_name[256];
  FILE *fp;
  int depth = 0;
  int *nx = grid_cell->nx[depth];
  int idx[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                   {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

  sprintf(file_name, "%s/grid.vtk", directory);
  fp = fopen(file_name, "w");

  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "Grid data\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fp, "POINTS %d float\n", nPoints);
  for (int ix = 0; ix < nx[0]; ix++) {
    for (int iy = 0; iy < nx[1]; iy++) {
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid_cell->grid[ix][iy][iz];

        for (int iDim = 0; iDim < 3; iDim++) {
          fprintf(fp, "%e ", igrid->x[iDim]);
        }
        fprintf(fp, "\n");

        if (igrid->subgrid != NULL) {
          WriteSubgridPoints(depth + 1, grid_cell, igrid->subgrid, fp);
        }
      }
    }
  }

  fprintf(fp, "CELLS %d %d\n", nCells, 9 * nCells);
  for (int ix = 0; ix < nx[0] - 1; ix++) {
    for (int iy = 0; iy < nx[1] - 1; iy++) {
      for (int iz = 0; iz < nx[2] - 1; iz++) {
        GRID_t *igrid = &grid_cell->grid[ix][iy][iz];

        fprintf(fp, "8 ");
        for (int i = 0; i < 8; i++) {
          fprintf(
              fp, "%d ",
              grid_cell->grid[ix + idx[i][0]][iy + idx[i][1]][iz + idx[i][2]]
                  .id);
        }
        fprintf(fp, "\n");

        if (igrid->subgrid != NULL) {
          WriteSubgridCells(depth + 1, grid_cell, igrid->subgrid, fp);
        }
      }
    }
  }

  fprintf (fp, "CELL_TYPES %d\n", nCells);
  for(int i = 0; i < nCells; i++) {
    fprintf (fp, "12\n");
  }

  fclose(fp);
}

int main(int argc, char **argv) {
  GRID_CELL_t grid_cell;

  ReadGrids(&grid_cell, argv[1]);
  WriteGrids(&grid_cell, argv[1]);

  return 0;
}