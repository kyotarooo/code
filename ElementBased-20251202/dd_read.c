#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dd_defs.h"
#include "dd_inclusion.h"
#include "dd_struct.h"

static void UnitVector(double* a) {
  double r = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    r += a[iDim] * a[iDim];
  }
  r = sqrt(r);
  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] /= r;
  }
}

static FILE* OpenFile(char* directory, char* file_name, char* mode) {
  char full_file_name[256];
  FILE* fp;

  sprintf(full_file_name, "%s/%s", directory, file_name);
  fp = fopen(full_file_name, mode);
  if (!strcmp(mode, "r") && fp == NULL) {
    fprintf(stderr, "Cannot open %s file\n", full_file_name);
    exit(1);
  }

  return fp;
}

static void ReadCondition(DD_t* dd, char* directory) {
  double g, v;
  double strain_rate;
  FILE* fp = OpenFile(directory, "condition.inp", "r");
  char time_integral[256];

  // Time stepping
  fscanf(fp, "%d", &dd->step.n);
  fscanf(fp, "%lf", &dd->step.umax);
  fscanf(fp, "%lf%lf", &dd->step.dtmin, &dd->step.dtmax);

  // Output
  fscanf(fp, "%d", &dd->output.interval);
  dd->output.id = 1;

  // Time integral method
  fscanf(fp, "%s", time_integral);

  if (!strcmp(time_integral, "trapezoid")) {
    dd->step.time_integral_method = TIME_INTEGRAL_TRAPEZOID;
  } else if (!strcmp(time_integral, "rk1")) {
    dd->step.time_integral_method = TIME_INTEGRAL_RK1;
  } else if (!strcmp(time_integral, "rk2")) {
    dd->step.time_integral_method = TIME_INTEGRAL_RK2;
  } else if (!strcmp(time_integral, "rk3")) {
    dd->step.time_integral_method = TIME_INTEGRAL_RK3;
  } else if (!strcmp(time_integral, "rk4")) {
    dd->step.time_integral_method = TIME_INTEGRAL_RK4;
  } else {
    fprintf(stderr, "Time integral method must be correctly selected.\n");
    exit(1);
  }

  // Trapesoidal time integral parameter
  fscanf(fp, "%lf", &dd->step.theta);

  // Initialize
  dd->step.t = 0.0;

  // Material constants
  fscanf(fp, "%lf%lf", &g, &v);
  dd->material.g = g;
  dd->material.v = v;
  dd->material.e = 2.0 * (1.0 + v) * g;
  fscanf(fp, "%lf", &dd->material.mobility);

  // Dislocation core radius
  fscanf(fp, "%lf", &dd->material.a_core);

  // Simulation volume
  dd->bc.size = (double*)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%lf", &dd->bc.size[iDim]);
  }

  // Periodic boundary condition
  dd->bc.periodic = (int*)malloc(3 * sizeof(int));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%d", &dd->bc.periodic[iDim]);
  }

  // Direct interaction distance
  fscanf(fp, "%lf", &dd->bc.r_interaction);

  // Applied stress (stress tensor must be symmetric)
  dd->bc.stress = (double**)malloc(3 * sizeof(double*));
  for (int iDim = 0; iDim < 3; iDim++) {
    dd->bc.stress[iDim] = (double*)malloc(3 * sizeof(double));
    for (int jDim = 0; jDim < 3; jDim++) {
      fscanf(fp, "%lf", &dd->bc.stress[iDim][jDim]);
    }
  }

  // Strain rate test
  dd->mechanical_behavior.index = (int*)malloc(2 * sizeof(int));
  dd->mechanical_behavior.strain_rate = (double**)malloc(3 * sizeof(double*));
  dd->mechanical_behavior.strain = (double**)malloc(3 * sizeof(double*));
  dd->mechanical_behavior.plastic_strain =
      (double**)malloc(3 * sizeof(double*));
  for (int iDim = 0; iDim < 3; iDim++) {
    dd->mechanical_behavior.strain_rate[iDim] =
        (double*)malloc(3 * sizeof(double));
    dd->mechanical_behavior.strain[iDim] = (double*)malloc(3 * sizeof(double));
    dd->mechanical_behavior.plastic_strain[iDim] =
        (double*)malloc(3 * sizeof(double));
    for (int jDim = 0; jDim < 3; jDim++) {
      dd->mechanical_behavior.strain_rate[iDim][jDim] = 0.0;
      dd->mechanical_behavior.strain[iDim][jDim] = 0.0;
      dd->mechanical_behavior.plastic_strain[iDim][jDim] = 0.0;
    }
  }
  for (int i = 0; i < 2; i++) {
    fscanf(fp, "%d", &dd->mechanical_behavior.index[i]);
  }
  fscanf(fp, "%lf", &strain_rate);
  if (dd->mechanical_behavior.index[0] != -1) {
    int i = dd->mechanical_behavior.index[0];
    int j = dd->mechanical_behavior.index[1];

    dd->mechanical_behavior.strain_rate[i][j] = strain_rate;
    dd->mechanical_behavior.strain_rate[j][i] = strain_rate;
  }

  // Standard element length
  dd->bc.element_length = (double*)malloc(2 * sizeof(double));
  for (int i = 0; i < 2; i++) {
    fscanf(fp, "%lf", &dd->bc.element_length[i]);
  }

  fclose(fp);
}

static void ReadDislocationCore(DD_t* dd, char* directory) {
  char file_name[256];
  FILE* fp;

  sprintf(file_name, "%s/core.inp", directory);
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    fprintf(stdout, "[LOG] No dislocation core data is found\n");
    fflush(stdout);
    dd->core.n = 0;
    dd->core.nDislocations = 0;
    return;
  }

  fprintf(stdout, "[LOG] Dislocation core data is found\n");
  fflush(stdout);

  fscanf(fp, "%d", &dd->core.n);
  fscanf(fp, "%d", &dd->core.nDislocations);

  dd->core.property = (PROPERTY_t*)malloc(dd->core.n * sizeof(PROPERTY_t));
  for (int iCore = 0; iCore < dd->core.n; iCore++) {
    PROPERTY_t* core = &dd->core.property[iCore];

    core->x = (double*)malloc(dd->core.nDislocations * sizeof(double));
    core->fmatrix = (double*)malloc(dd->core.nDislocations * sizeof(double));
    core->finclusion = (double*)malloc(dd->core.nDislocations * sizeof(double));
    core->bedge = (double*)malloc(dd->core.nDislocations * sizeof(double));
    core->bscrew = (double*)malloc(dd->core.nDislocations * sizeof(double));
    for (int i = 0; i < dd->core.nDislocations; i++) {
      fscanf(fp, "%lf%lf%lf%lf%lf", &core->x[i], &core->fmatrix[i],
             &core->finclusion[i], &core->bedge[i], &core->bscrew[i]);
    }
  }
  fclose(fp);
}

static void ReadElements(DD_t* dd, char* directory) {
  FILE* fp = OpenFile(directory, "elem.inp", "r");

  dd->node = (NODE_t*)malloc(NMAX_NODES * sizeof(NODE_t));
  dd->element = (ELEMENT_t*)malloc(NMAX_ELEMENTS * sizeof(ELEMENT_t));
  for (int iNode = 0; iNode < NMAX_NODES; iNode++) {
    NODE_t* node = &dd->node[iNode];

    node->elementID = (int*)malloc(NMAX_ELEMENTS_AROUND_NODE * sizeof(int));
    node->x = (double*)malloc(3 * sizeof(double));
    node->x0 = (double*)malloc(3 * sizeof(double));
    node->v = (double*)malloc(3 * sizeof(double));
    node->v0 = (double*)malloc(3 * sizeof(double));
    node->v1 = (double*)malloc(3 * sizeof(double));
    node->v2 = (double*)malloc(3 * sizeof(double));
    node->v3 = (double*)malloc(3 * sizeof(double));
    node->direction = (double*)malloc(3 * sizeof(double));
  }
  for (int iElement = 0; iElement < NMAX_ELEMENTS; iElement++) {
    ELEMENT_t* element = &dd->element[iElement];

    element->nodeID = (int*)malloc(2 * sizeof(int));
    element->burgers = (double*)malloc(3 * sizeof(double));
    element->slip = (double*)malloc(3 * sizeof(double));
    element->tangent = (double*)malloc(3 * sizeof(double));
    element->b = (double*)malloc(6 * sizeof(double));
    element->K = (double**)malloc(6 * sizeof(double*));
    for (int i = 0; i < 6; i++) {
      element->K[i] = (double*)malloc(6 * sizeof(double));
    }
    element->integral =
        (INTEGRAL_t*)malloc(N_DD_INTEGRALS * sizeof(INTEGRAL_t));
    if (dd->core.nDislocations != 0) {
      for (int iIntegral = 0; iIntegral < N_DD_INTEGRALS; iIntegral++) {
        INTEGRAL_t* integral = &element->integral[iIntegral];

        integral->core_x =
            (double*)malloc(dd->core.nDislocations * sizeof(double));
      }
    }
    // Default value
    element->coreID = 0;
  }

  fscanf(fp, "%d", &dd->nNodes);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t* node = &dd->node[iNode];

    fscanf(fp, "%d", &node->type);
    for (int iDim = 0; iDim < 3; iDim++) {
      fscanf(fp, "%lf", &node->x[iDim]);
      node->direction[iDim] = 0.0;  // Initial setup
    }
  }

  fscanf(fp, "%d", &dd->nElements);
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t* element = &dd->element[iElement];

    for (int iNode = 0; iNode < 2; iNode++) {
      fscanf(fp, "%d", &element->nodeID[iNode]);
      if (element->nodeID[iNode] >= dd->nNodes) {
        fprintf(stdout, "Node ID is larger than the maximum\n");
        fprintf(stdout, "Node ID: %d\n", element->nodeID[iNode]);
        fprintf(stdout, "Maximum node ID: %d\n", dd->nNodes - 1);
        exit(1);
      }
    }
    for (int iDim = 0; iDim < 3; iDim++) {
      fscanf(fp, "%lf", &element->burgers[iDim]);
    }
    for (int iDim = 0; iDim < 3; iDim++) {
      fscanf(fp, "%lf", &element->slip[iDim]);
    }
    UnitVector(element->slip);
    // If you have dislocation cores, additional data must be read.
    if (dd->core.nDislocations != 0) {
      int coreID;

      fscanf(fp, "%d", &coreID);
      element->coreID = coreID;
      // Substitute the initial dislocation core structure
      for (int iIntegral = 0; iIntegral < N_DD_INTEGRALS; iIntegral++) {
        INTEGRAL_t* integral = &element->integral[iIntegral];

        for (int i = 0; i < dd->core.nDislocations; i++) {
          integral->core_x[i] = dd->core.property[coreID].x[i];
        }
      }
    }
  }

  fclose(fp);
}

static void ReadGridArrangement(DD_t* dd, char* directory) {
  INITIAL_STRESS_t* initial_stress = &dd->initial_stress;
  FILE* fp;
  char file_name[256];

  initial_stress->nx = (int*)malloc(3 * sizeof(int));
  initial_stress->size = (double*)malloc(3 * sizeof(double));

  // Initialize
  for (int iDim = 0; iDim < 3; iDim++) {
    initial_stress->nx[iDim] = 1;
    initial_stress->size[iDim] = dd->bc.size[iDim];
  }

  sprintf(file_name, "%s/grid_arrangement.inp", directory);
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    fprintf(stdout, "[LOG] No grid arrangement data is found\n");
    fflush(stdout);
    initial_stress->have_grids = 0;
    return;
  }

  fprintf(stdout, "[LOG] Grid arrangement data is found\n");
  fflush(stdout);

  initial_stress->have_grids = 1;

  for(int iDim = 0; iDim < 3; iDim++) {
    fscanf (fp, "%lf", &initial_stress->size[iDim]);
  }

  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%d", &initial_stress->nx[iDim]);
  }

  // Scaling factor to control eigen strain value
  fscanf (fp, "%lf", &initial_stress->r);

  fclose(fp);
}

static void ReadSubgrid(int depth, GRID_CELL_t* grid_cell, GRID_t* grid,
                        double r, FILE* fp) {
  int* nx = grid_cell->nx[depth];

  grid->subgrid = (GRID_t***)malloc(nx[0] * sizeof(GRID_t**));
  for (int ix = 0; ix < nx[0]; ix++) {
    grid->subgrid[ix] = (GRID_t**)malloc(nx[1] * sizeof(GRID_t*));
    for (int iy = 0; iy < nx[1]; iy++) {
      grid->subgrid[ix][iy] = (GRID_t*)malloc(nx[2] * sizeof(GRID_t));
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t* subgrid = &grid->subgrid[ix][iy][iz];
        int have_subgrid;

        subgrid->stress = (double**)malloc(3 * sizeof(double*));
        subgrid->stress[0] = (double*)malloc(9 * sizeof(double));
        for (int iDim = 1; iDim < 3; iDim++) {
          subgrid->stress[iDim] = subgrid->stress[iDim - 1] + 3;
        }
        fread(subgrid->stress[0], sizeof(double), 9, fp);

        for(int iDim = 0; iDim < 3; iDim++) {
          for(int jDim = 0; jDim < 3; jDim++) {
            subgrid->stress[iDim][jDim] *= r;
          }
        }

        fread(&have_subgrid, sizeof(int), 1, fp);
        if (have_subgrid) {
          ReadSubgrid(depth + 1, grid_cell, subgrid, r, fp);
        } else {
          subgrid->subgrid = NULL;
        }
      }
    }
  }
}

static void ReadGridStress(DD_t* dd, char* directory) {
  char file_name[256];
  FILE* fp;
  int depth = 0;
  INITIAL_STRESS_t* initial_stress = &dd->initial_stress;
  GRID_CELL_t* grid_cell = &initial_stress->grid_cell;
  double r = initial_stress->r;
  int* nx;

  sprintf(file_name, "%s/grid_stress.inp", directory);
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    fprintf(stdout, "[LOG] No grid stress data is found\n");
    fflush(stdout);
    return;
  }

  fprintf(stdout, "[LOG] Grid stress data is found\n");
  fflush(stdout);

  fread(&grid_cell->max_depth, sizeof(int), 1, fp);
  fprintf (stdout, "[LOG] Grid depth: %d\n", grid_cell->max_depth);
  grid_cell->nx = (int**)malloc(grid_cell->max_depth * sizeof(int*));
  for (int i = 0; i < grid_cell->max_depth; i++) {
    fprintf (stdout, "[LOG] Number of grid points (%d): ", i + 1);
    grid_cell->nx[i] = (int*)malloc(3 * sizeof(int));
    fread(grid_cell->nx[i], sizeof(int), 3, fp);
    for(int iDim = 0; iDim < 3; iDim++) {
      fprintf (stdout, "%d ", grid_cell->nx[i][iDim]);
    }
    fprintf (stdout, "\n");
  }
  fflush(stdout);
  nx = grid_cell->nx[depth];

  grid_cell->size = (double*)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    grid_cell->size[iDim] =
        initial_stress->size[iDim] / ((double)nx[iDim] - 1.0);
  }

  grid_cell->grid = (GRID_t***)malloc(nx[0] * sizeof(GRID_t**));
  for (int ix = 0; ix < nx[0]; ix++) {
    grid_cell->grid[ix] = (GRID_t**)malloc(nx[1] * sizeof(GRID_t*));
    for (int iy = 0; iy < nx[1]; iy++) {
      grid_cell->grid[ix][iy] = (GRID_t*)malloc(nx[2] * sizeof(GRID_t));
      for (int iz = 0; iz < nx[2]; iz++) {
        GRID_t* grid = &grid_cell->grid[ix][iy][iz];
        int have_subgrid;

        grid->stress = (double**)malloc(3 * sizeof(double*));
        grid->stress[0] = (double*)malloc(9 * sizeof(double));
        for (int iDim = 1; iDim < 3; iDim++) {
          grid->stress[iDim] = grid->stress[iDim - 1] + 3;
        }
        fread(grid->stress[0], sizeof(double), 9, fp);

        for(int iDim = 0; iDim < 3; iDim++) {
          for(int jDim = 0; jDim < 3; jDim++) {
            grid->stress[iDim][jDim] *= r;
          }
        }
        
        fread(&have_subgrid, sizeof(int), 1, fp);
        if (have_subgrid) {
          ReadSubgrid(depth + 1, grid_cell, grid, r, fp);
        } else {
          grid->subgrid = NULL;
        }
      }
    }
  }

  fclose(fp);
}

static void ReadInclusions(DD_t* dd, char* directory) {
  char file_name[256];
  FILE* fp;

  sprintf(file_name, "%s/inclusion.inp", directory);
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    fprintf(stdout, "[LOG] No inclusion data is found\n");
    fflush(stdout);
    dd->nInclusions = 0;
    return;
  }

  fprintf(stdout, "[LOG] Inclusion data is found\n");
  fflush(stdout);

  fscanf(fp, "%d", &dd->nInclusions);
  dd->inclusion = (INCLUSION_t*)malloc(dd->nInclusions * sizeof(INCLUSION_t));
  for (int i = 0; i < dd->nInclusions; i++) {
    INCLUSION_t* inclusion = &dd->inclusion[i];

    fscanf(fp, "%d", &inclusion->type);

    inclusion->x = (double*)malloc(3 * sizeof(double));
    for (int iDim = 0; iDim < 3; iDim++) {
      fscanf(fp, "%lf", &inclusion->x[iDim]);
    }

    switch (inclusion->type) {
      case SPHERICAL_INCLUSION:
        inclusion->v = (double*)malloc(1 * sizeof(double));
        fscanf(fp, "%lf", &inclusion->v[0]);  // Inclusion radius
        break;
      case CYLINDERICAL_INCLUSION:
        inclusion->v = (double*)malloc(5 * sizeof(double));
        // Radius, height, and major axis(x y z)
        for (int i = 0; i < 5; i++) {
          fscanf(fp, "%lf", &inclusion->v[i]);
        }
        DD_SetUpCylindericalInclusionCoordinateSystem(inclusion);
        break;
    }
  }
  fclose(fp);
}

void DD_ReadInput(DD_t* dd, char* directory) {
  fprintf(stdout, "[LOG] Reading DD inputs\n");
  fflush(stdout);
  ReadCondition(dd, directory);
  ReadDislocationCore(dd, directory);
  ReadElements(dd, directory);
  ReadGridArrangement(dd, directory);
  if (dd->initial_stress.have_grids) {
    ReadGridStress(dd, directory);
  }
  ReadInclusions(dd, directory);
}
