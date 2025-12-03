#ifndef __FEM_STRUCT_H_
#define __FEM_STRUCT_H_

#ifndef _NOFEM_
#ifdef _PARDISO_
#include <mkl.h>
#else
#include "dmumps_c.h"
#endif
#endif

typedef struct __fem_node_t_ {
  int iTmp;
  int *bc;
  double *x;
  double *u;
  double *bc_u;
  double *f;
  double *bc_f;
  double *stress;
} FEMNODE_t;

typedef struct __fem_element_t_ {
  int *nodeID;
  int nSurfaces;
  int *surface;
  double *f, *f_grid;
  double *xmin, *xmax;
  double **surface_nx;
  double **surface_x;
  double *x;
} FEMELEMENT_t;

typedef struct __fem_material_t_ {
  double g, v, e;
} FEMMATERIAL_t;

typedef struct __fem_equation_t_ {
  int n;
  int *nNodes, **nodeID;
  int *start_index;
  double *a, *x, *b;
  double scale_u, scale_f;
  
#ifndef _NOFEM_
#ifdef _PARDISO_
  int *bc;
  double *b0;
  int *nDofs, **dofID;
  MKL_INT *MKL_system_index;
  MKL_INT *MKL_start_index;
#else
  int *row, *col;
  double *bc;
  double *a0;
  DMUMPS_STRUC_C id;
#endif
#endif
} FEMEQUATION_t;

typedef struct __fem_t_ {
  int nNodes;
  FEMNODE_t *node;
  int nElements;
  FEMELEMENT_t *element;
  FEMMATERIAL_t material;
  FEMEQUATION_t equation;

  double *xmin, *xmax;

  int nSurface_elements;
  int *surface_element;
  int **surface_surface;

  double volume;
} FEM_t;

#endif
