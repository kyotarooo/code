#ifndef __DD_STRUCT_H_
#define __DD_STRUCT_H_

#include <stdio.h>

typedef struct __node_t_ {
  int type;                       // Node type
  int nElements;                  // Number of elements around the node
  int* elementID;                 // Element IDs
  double *x, *x0;                 // Position vector
  double *v, *v0, *v1, *v2, *v3;  // Velocity vectors for Runge-Kutta integral
                                  // Working vectors for Trapezoidal integral
  double* direction;  // If this node is "constraint node", the node can move
                      // only in this direction
} NODE_t;

typedef struct __integral_t_ {
  double* core_x;  // 1-D coordinate of virtual dislocation structure
} INTEGRAL_t;

typedef struct __element_t_ {
  int iTmp;
  int* nodeID;           // Node IDs consisting of the element
  double* burgers;       // Burgers vector
  double* slip;          // Slip plane unit normal vector
  double* tangent;       // Unit tangent vector
  double length;         // Element length
  double **K, *b;        // Coefficient and RHS of global equation
  int coreID;            // Dislocation core type
  INTEGRAL_t* integral;  // Dislocation core structure
} ELEMENT_t;

typedef struct __material_t_ {
  double e, g, v;   // Young's modulus, Poisson's ratio, shear modulus
  double mobility;  // Mobility
  double a_core;    // Dislocation core radius
} MATERIAL_t;

typedef struct __step_t_ {
  int n;                // Number of time steps
  double t, dt;         // Time, time increment
  double umax;          // Maximum displacement per time step
  double dtmin, dtmax;  // Minimum and maximum time increment

  int time_integral_method;
  double theta;  // Control the time integral scheme (explicit, trapezoidal or
                 // implicit)
} STEP_t;

typedef struct __output_t_ {
  int interval;  // Output interval
  int id;        // Output ID
} OUTPUT_t;

typedef struct __bc_t_ {
  double* size;   // Size of simulation volume in each direction
  int* periodic;  // Periodic boundary condition (1 or 0)
  double*
      element_length;    // Standard element length to control element divisions
  double** stress;       // Applied stresss tensor
  double r_interaction;  // Critical distance for direct interaction between
                         // elements
} BC_t;

typedef struct __mechanical_behavior_t_ {
  int* index;               // Index of strain/stress tensor
  double** strain_rate;     // Strain rate
  double** strain;          // Accumulated strain
  double** plastic_strain;  // Accumurated plastic strain
} MECHANICAL_BEHAVIOR_t;

typedef struct __subgrid_stress_t_ {
  int n;              // Number of subgrids
  int* nx;            // Number of subgrids in each axis
  double* size;       // Size of each subgrid
  int* gridID;        // Parent grid ID
  double**** stress;  // Stress tensr at subgrid points
} SUBGRID_t;

typedef struct __grid_t_ {
  double** stress;
  struct __grid_t_*** subgrid;
} GRID_t;

typedef struct __grid_cell_t_ {
  int max_depth;
  int** nx;
  double* size;
  GRID_t*** grid;
} GRID_CELL_t;

typedef struct __initial_stress_t_ {
  int have_grids;         // If we have grids
  int* nx;                // Number of grid-cells in each axis
  double* size;           // Size of each grid-cell in each direction
  GRID_CELL_t grid_cell;  // Grid data
  double r;               // Scaling factor to control eigen strain value
} INITIAL_STRESS_t;

typedef struct __inclusion_t_ {
  int type;   // Type
  double* x;  // Position
  double* v;  // Parameters
  double* vector;
  double** tensor;
} INCLUSION_t;

typedef struct __property_t_ {
  double* x;                     // Templete dislocation core structure
  double *fmatrix, *finclusion;  // Lattice restoring force
  double *bedge, *bscrew;        // burgers vectors
} PROPERTY_t;

typedef struct __core_t_ {
  int n;                 // Number of dislocation core types
  int nDislocations;     // Number of fractional dislocations
  PROPERTY_t* property;  // Dislocation core property data
} CORE_t;

typedef struct __file_t_ {
  FILE *direct_interaction;
} FILE_t;

typedef struct __dd_t_ {
  int nNodes;
  NODE_t* node;
  int nElements;
  ELEMENT_t* element;
  MATERIAL_t material;
  STEP_t step;
  OUTPUT_t output;
  BC_t bc;
  MECHANICAL_BEHAVIOR_t mechanical_behavior;
  INITIAL_STRESS_t initial_stress;
  int nInclusions;
  INCLUSION_t* inclusion;
  CORE_t core;
  FILE_t file;
} DD_t;

#endif
