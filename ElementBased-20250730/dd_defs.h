#ifndef __DD_DEFS_H_
#define __DD_DEFS_H_

// Log inverval
#define LOG_INTERVAL 100

// Maximum number of nodes
#define NMAX_NODES 1000000
// Maximum number of elements
#define NMAX_ELEMENTS 1000000
// Maximum number of element pairs
#define NMAX_ELEMENT_PAIRS 10

// Maximum number of elements around a node
#define NMAX_ELEMENTS_AROUND_NODE 10

// Mobile node
#define MOBILE_NODE 0
// Immobile node
#define IMMOBILE_NODE 1
// Constraint node
#define CONSTRAINT_NODE 2
// Node on free surface
#define FREE_SURFACE_NODE 3

// NULL flag
#define NO_NODE -1
#define NO_ELEMENT -1

// Number of integrals for velocity calculation
#define N_DD_INTEGRALS 4

// Upper limit of cosine angle between elements
#define SMOOTH_COSINE_ANGLE 0.7

// lower limit of cosine angle between element and free surface
#define FREE_SURFACE_COSINE_ANGLE 0.7

// Trapezoidal time integration
#define NMAX_TRAPEZOIDAL_ITERATIONS 10
#define TRAPEZOIDAL_EPS 1.0e-04
#define N_TRAPEZOIDAL_SUBDIVISIONS 2

// Direct interactions
#define NO_INTERACTION -1
#define ANNIHILATION 0
#define JUNCTION 1
#define DIRECT_INTERACTION_COSINE_ANGLE 0.8

// Interval for update of correction field
#define CORRECTION_FIELD_UPDATE_INTERVAL 100

// Subgrid
#define NO_SUBGRID -1

// Restart
#define READ_RESTART_DATA 0
#define WRITE_RESTART_DATA 1

// Path to gnuplot
#define GNUPLOT_PATH "/opt/homebrew/bin/gnuplot"

#endif