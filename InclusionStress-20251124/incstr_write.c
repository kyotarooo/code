#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "incstr_grid.h"
#include "incstr_struct.h"

static void WriteSubgridStress(int depth, GRID_CELL_t *grid_cell,
                               GRID_t ***grid, FILE *fp) {
  int *nx = grid_cell->nx[depth];
  double stress_vector[9];
  int index[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

  for (int ix = 0; ix < nx[0]; ix++) {
    for (int iy = 0; iy < nx[1]; iy++) {
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t *igrid = &grid[ix][iy][iz];
        int have_subgrid = 0;

        for(int i = 0; i < 9; i++) {
          stress_vector[i] = igrid->stress[index[i]];
        }
        fwrite(stress_vector, sizeof(double), 9, fp);

        if (igrid->subgrid != NULL) have_subgrid = 1;
        fwrite(&have_subgrid, sizeof(int), 1, fp);

        if(have_subgrid) {
          WriteSubgridStress(depth + 1, grid_cell, igrid->subgrid, fp);
        }
      }
    }
  }
}

void INCSTR_WriteStress(GRID_CELL_t *grid_cell, char *directory) {
  char file_name[256];
  FILE *fp;
  double stress_vector[9];
  int index[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

  sprintf(file_name, "%s/grid_stress.inp", directory);
  fp = fopen(file_name, "w");

  fwrite(&grid_cell->max_depth, sizeof(int), 1, fp);
  for (int i = 0; i < grid_cell->max_depth; i++) {
    fwrite(grid_cell->nx[i], sizeof(int), 3, fp);
  }

  for (int ix = 0; ix < grid_cell->nx[0][0]; ix++) {
    for (int iy = 0; iy < grid_cell->nx[0][1]; iy++) {
      for (int iz = 0; iz < grid_cell->nx[0][2]; iz++) {
        GRID_t *igrid = &grid_cell->grid[ix][iy][iz];
        int have_subgrid = 0;

        for (int i = 0; i < 9; i++) {
          stress_vector[i] = igrid->stress[index[i]];
        }
        fwrite(stress_vector, sizeof(double), 9, fp);

        if (igrid->subgrid != NULL) have_subgrid = 1;
        fwrite(&have_subgrid, sizeof(int), 1, fp);

        if (have_subgrid) {
          WriteSubgridStress(1, grid_cell, igrid->subgrid, fp);
        }
      }
    }
  }

  fclose(fp);
}