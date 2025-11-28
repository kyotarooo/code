#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "incstr_defs.h"
#include "incstr_math.h"
#include "incstr_struct.h"

static FILE* OpenFile(char* file_name, char* directory) {
  char full_file_name[256];
  FILE* fp;

  sprintf(full_file_name, "%s/%s", directory, file_name);
  fp = fopen(full_file_name, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot open %s file\n", full_file_name);
    exit(1);
  }
  return fp;
}

static void TinyNoise(double r, double* x) {
  double ratio = 1.0e-06;

  for (int iDim = 0; iDim < 3; iDim++) {
    double noise = ratio * r * ((double)rand() / (double)RAND_MAX - 0.5);

    x[iDim] += noise;
  }
}

void INCSTR_ReadInclusions(INCLUSION_t* inclusion, char* directory) {
  FILE* fp;

  fp = OpenFile("inclusion.inp", directory);

  fscanf(fp, "%d", &inclusion->n);

  inclusion->shape = (int*)malloc(inclusion->n * sizeof(int));
  inclusion->eigen_strain = (double*)malloc(inclusion->n * sizeof(double));
  inclusion->size = (double**)malloc(inclusion->n * sizeof(double*));
  inclusion->x = (double**)malloc(inclusion->n * sizeof(double*));
  inclusion->orientation = (double***)malloc(inclusion->n * sizeof(double**));
  inclusion->elastic_dipole = (double**)malloc(inclusion->n * sizeof(double*));

  for (int i = 0; i < inclusion->n; i++) {
    char shape[256];

    fscanf(fp, "%s", shape);
    fscanf(fp, "%lf", &inclusion->eigen_strain[i]);

    if (!strcmp(shape, "sphere")) {
      inclusion->shape[i] = SPHERE;
      inclusion->size[i] = (double*)malloc(1 * sizeof(double));
      fscanf(fp, "%lf", &inclusion->size[i][0]);  // Radius
    } else if (!strcmp(shape, "cylinder")) {
      inclusion->shape[i] = CYLINDER;
      inclusion->size[i] = (double*)malloc(3 * sizeof(double));
      fscanf(fp, "%lf", &inclusion->size[i][0]);  // Radius
      fscanf(fp, "%lf", &inclusion->size[i][1]);  // Length of major axis
      fscanf(fp, "%lf", &inclusion->size[i][2]);  // 1st nnbr distance
    } else if (!strcmp(shape, "axial_cylinder")) {
      inclusion->shape[i] = AXIAL_CYLINDER;
      inclusion->size[i] = (double*)malloc(3 * sizeof(double));
      fscanf(fp, "%lf", &inclusion->size[i][0]);  // Radius
      fscanf(fp, "%lf", &inclusion->size[i][1]);  // Length of major axis
      fscanf(fp, "%lf", &inclusion->size[i][2]);  // 1st nnbr distance
    } else if (!strcmp(shape, "cuboid")) {
      inclusion->shape[i] = CUBOID;
      inclusion->size[i] = (double*)malloc(3 * sizeof(double));
      fscanf(fp, "%lf", &inclusion->size[i][0]);  // Edge size
      fscanf(fp, "%lf", &inclusion->size[i][1]);  // Edge size
      fscanf(fp, "%lf", &inclusion->size[i][2]);  // Edge size
    } else if (!strcmp(shape, "truncated_sphere")) {
      inclusion->shape[i] = TRUNCATED_SPHERE;
      inclusion->size[i] = (double*)malloc(3 * sizeof(double));
      fscanf(fp, "%lf", &inclusion->size[i][0]);  // Radius of sphere
      fscanf(fp, "%lf", &inclusion->size[i][1]);  // Bottom position
      fscanf(fp, "%lf", &inclusion->size[i][2]);  // Top position
    } else if (!strcmp(shape, "elastic_dipole")) {
      inclusion->shape[i] = ELASTIC_DIPOLE;
      inclusion->size[i] = (double*)malloc(1 * sizeof(double));
      fscanf(fp, "%lf",
             &inclusion->size[i][0]);  // Core radius for non-singular
    } else if (!strcmp(shape, "dilatation_center")) {
      inclusion->shape[i] = DILATATION_CENTER;
      inclusion->size[i] = (double*)malloc(1 * sizeof(double));
      fscanf(fp, "%lf", &inclusion->size[i][0]);  // Lattice constant
    }

    inclusion->x[i] = (double*)malloc(3 * sizeof(double));
    for (int iDim = 0; iDim < 3; iDim++) {
      fscanf(fp, "%lf", &inclusion->x[i][iDim]);
    }

    if (!strcmp(shape, "cylinder") || !strcmp(shape, "axial_cylinder")) {
      inclusion->orientation[i] = (double**)malloc(1 * sizeof(double*));
      inclusion->orientation[i][0] = (double*)malloc(3 * sizeof(double));
      for (int iDim = 0; iDim < 3; iDim++) {
        fscanf(fp, "%lf", &inclusion->orientation[i][0][iDim]);
      }
      INCSTR_NormalizeVector(inclusion->orientation[i][0]);
      TinyNoise(inclusion->size[i][0], inclusion->x[i]);
    } else if (!strcmp(shape, "cuboid")) {
      inclusion->orientation[i] = (double**)malloc(3 * sizeof(double));
      for (int iDim = 0; iDim < 3; iDim++) {
        inclusion->orientation[i][iDim] = (double*)malloc(3 * sizeof(double));
        for (int jDim = 0; jDim < 3; jDim++) {
          fscanf(fp, "%lf", &inclusion->orientation[i][iDim][jDim]);
        }
      }
      INCSTR_NormalizeTensor(inclusion->orientation[i]);
      TinyNoise(inclusion->size[i][0], inclusion->x[i]);
    } else if (!strcmp(shape, "truncated_sphere")) {
      inclusion->orientation[i] = (double**)malloc(1 * sizeof(double*));
      inclusion->orientation[i][0] = (double*)malloc(3 * sizeof(double));
      for (int iDim = 0; iDim < 3; iDim++) {
        fscanf(fp, "%lf", &inclusion->orientation[i][0][iDim]);
      }
      INCSTR_NormalizeVector(inclusion->orientation[i][0]);
      // To avoid any singularly in the analytical solution,
      // the position of inclusion is slightly moved.
      TinyNoise(inclusion->size[i][0], inclusion->x[i]);
    } else if (!strcmp(shape, "elastic_dipole")) {
      int idx = 0;

      inclusion->elastic_dipole[i] = (double*)malloc(9 * sizeof(double));
      for (int iDim = 0; iDim < 3; iDim++) {
        for (int jDim = 0; jDim < 3; jDim++) {
          fscanf(fp, "%lf", &inclusion->elastic_dipole[i][idx]);
          idx = idx + 1;
        }
      }
      TinyNoise(1.0e-10, inclusion->x[i]);
    } else if (!strcmp(shape, "dilatation_center")) {
      TinyNoise(1.0e-10, inclusion->x[i]);
    }
  }

  fclose(fp);
}

void INCSTR_ReadGrids(GRID_CELL_t* grid_cell, MATERIAL_t* material,
                      char* directory) {
  FILE* fp = OpenFile("grid.inp", directory);
  double* size;

  fscanf(fp, "%d", &grid_cell->max_depth);
  grid_cell->nx = (int**)malloc(grid_cell->max_depth * sizeof(int*));
  grid_cell->size = (double**)malloc(grid_cell->max_depth * sizeof(double*));
  size = material->size;
  for (int i = 0; i < grid_cell->max_depth; i++) {
    grid_cell->nx[i] = (int*)malloc(3 * sizeof(int));
    grid_cell->size[i] = (double*)malloc(3 * sizeof(double));

    for (int iDim = 0; iDim < 3; iDim++) {
      fscanf(fp, "%d", &grid_cell->nx[i][iDim]);  // Number of grids
      grid_cell->nx[i][iDim] += 1;                // Number of grid points
      grid_cell->size[i][iDim] =
          size[iDim] / (double)(grid_cell->nx[i][iDim] - 1);
    }
    size = grid_cell->size[i];
  }

  fclose(fp);
}

void INCSTR_ReadMaterial(MATERIAL_t* material, char* directory) {
  FILE* fp = OpenFile("material.inp", directory);

  material->orientation = (double**)malloc(3 * sizeof(double*));
  for (int iDim = 0; iDim < 3; iDim++) {
    material->orientation[iDim] = (double*)malloc(3 * sizeof(double));
    for (int jDim = 0; jDim < 3; jDim++) {
      fscanf(fp, "%lf", &material->orientation[iDim][jDim]);
    }
  }
  INCSTR_NormalizeTensor(material->orientation);

  material->periodic = (int*)malloc(3 * sizeof(int));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%d", &material->periodic[iDim]);
  }

  material->size = (double*)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%lf", &material->size[iDim]);
  }

  fscanf(fp, "%lf%lf", &material->g, &material->v);

  fclose(fp);
}