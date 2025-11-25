#ifndef __INCSTR_STRUCT_H_
#define __INCSTR_STRUCT_H_

typedef struct __inclusion_t_ {
  int n;                    // Number of inclusions
  int *shape;               // Inclusion shapes
  double *eigen_strain;     // Eigen strain
  double **size;            // Parameters associated with each inclusion type
  double **x;               // Center position
  double ***orientation;    // Crystal orientation of inclusion
  double **elastic_dipole;  // Elastic dipole tensor
} INCLUSION_t;

typedef struct __grid_t_ {
  int inside;
  double *x;
  double *stress;
  struct __grid_t_ ***subgrid;
} GRID_t;

typedef struct __grid_cell_t_ {
  int max_depth;
  int **nx;
  double **size;
  GRID_t ***grid;
} GRID_CELL_t;

typedef struct __material_t_ {
  double **orientation;  // Crystal orientation of each axis
  int *periodic;         // Periodic boundary condition
  double *size;          // Simulation volume size
  double g, v;           // Elastic shear modulus and Poisson's ratio
  int nImages;      // Number of image cells for the application of periodic
                    // boundary condition
  double **offset;  // Offset position of the image cells
} MATERIAL_t;

#endif