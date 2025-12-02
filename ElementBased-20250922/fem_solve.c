#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fem_struct.h"

#ifndef _NO_FEM_
#ifdef _PARDISO_
#include <mkl.h>

#define PARDISO pardiso_

static void *MKL_pt[64];
static int MKL_mtype = -2; /* Real and unsymmetric */
static int MKL_nrhs = 1;   /* number of RHS */
static int MKL_maxfct = 1; /* maximum number of numerical factrizations */
static int MKL_mnum = 1;   /* which factrization to use */
static int MKL_iparm[64];
static int MKL_msglvl = 0;

static void MakeUnknown(FEM_t *fem) {
  double *x = fem->equation.x;
  double *b = fem->equation.b0;

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];
    for (int iDim = 0; iDim < 3; iDim++) {
      x[3 * iNode + iDim] = node->bc_u[iDim];
      b[3 * iNode + iDim] = 0.0;
    }
  }
}

static void SetDirichletBCForMatrix(int n, FEMEQUATION_t *equation) {
  int *bc = equation->bc;
  double *a = equation->a;
  double *x = equation->x;
  double *b = equation->b0;
  int *start_index = equation->start_index;
  int *nDofs = equation->nDofs;
  int **dofID = equation->dofID;

  for (int iDof = 0; iDof < n; iDof++) {
    if (bc[iDof]) {
      int start_indexi = start_index[iDof];

      for (int jDof = 0; jDof < nDofs[iDof]; jDof++) {
        int jDofID = dofID[iDof][jDof];

        b[jDofID] -= a[start_indexi + jDof] * x[iDof];
      }
    }
  }

  for (int iDof = 0; iDof < n; iDof++) {
    int start_indexi = start_index[iDof];
    for (int jDof = 0; jDof < nDofs[iDof]; jDof++) {
      int jDofID = dofID[iDof][jDof];
      if (bc[jDofID] && iDof < jDofID) {
        b[iDof] -= a[start_indexi + jDof] * x[jDofID];
      }
    }
  }

  for (int iDof = 0; iDof < n; iDof++) {
    int start_indexi = start_index[iDof];
    for (int jDof = 0; jDof < nDofs[iDof]; jDof++) {
      int jDofID = dofID[iDof][jDof];
      if (bc[jDofID]) {
        a[start_indexi + jDof] = 0.0;
      }
    }
  }

  for (int iDof = 0; iDof < n; iDof++) {
    if (bc[iDof]) {
      int start_indexi = start_index[iDof];

      a[start_indexi + 0] = 1.0;
      for (int jDof = 1; jDof < nDofs[iDof]; jDof++) {
        a[start_indexi + jDof] = 0.0;
      }
    }
  }
}

void FEM_Factorization(FEM_t *fem) {
  int nNodes = fem->nNodes;
  int n = 3 * nNodes;
  FEMEQUATION_t *equation = &fem->equation;
  int phase, error, idum;
  double ddum;

  if (fem->nElements == 0) return;

  fprintf(stdout, "[LOG] Compute the inverse of the stiffness matrix\n");
  fflush(stdout);

  MakeUnknown(fem);
  SetDirichletBCForMatrix(n, equation);

  MKL_iparm[0] = 0;
  for (int i = 0; i < 64; i++) {
    MKL_pt[i] = 0;
  }

  phase = 11;
  error = 0;
  idum = 0;
  ddum = 0.0;
  PARDISO(MKL_pt, &MKL_maxfct, &MKL_mnum, &MKL_mtype, &phase, &n, equation->a,
          equation->MKL_start_index, equation->MKL_system_index, &idum,
          &MKL_nrhs, MKL_iparm, &MKL_msglvl, &ddum, &ddum, &error);

  phase = 22;
  error = 0;
  idum = 0;
  ddum = 0.0;
  PARDISO(MKL_pt, &MKL_maxfct, &MKL_mnum, &MKL_mtype, &phase, &n, equation->a,
          equation->MKL_start_index, equation->MKL_system_index, &idum,
          &MKL_nrhs, MKL_iparm, &MKL_msglvl, &ddum, &ddum, &error);
}

static void MakeRHS(FEM_t *fem) {
  double *b = fem->equation.b;
  double *x = fem->equation.x;
  double *b0 = fem->equation.b0;
  double scale_u = fem->equation.scale_u;
  double scale_f = fem->equation.scale_f;

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    for (int iDim = 0; iDim < 3; iDim++) {
      b[3 * iNode + iDim] = scale_u * b0[3 * iNode + iDim] + node->f[iDim] +
                            scale_f * node->bc_f[iDim];
      x[3 * iNode + iDim] = 0.0;
    }
  }
}

void FEM_ComputeSolution(FEM_t *fem) {
  int nNodes = fem->nNodes;
  int n = 3 * nNodes;
  FEMEQUATION_t *equation = &fem->equation;
  int phase, error, idum;
  double ddum;
  double scale_u = equation->scale_u;

  if (fem->nElements == 0) return;

  MakeRHS(fem);

  phase = 33;
  error = 0;
  idum = 0;
  ddum = 0.0;
  PARDISO(MKL_pt, &MKL_maxfct, &MKL_mnum, &MKL_mtype, &phase, &n, equation->a,
          equation->MKL_start_index, equation->MKL_system_index, &idum,
          &MKL_nrhs, MKL_iparm, &MKL_msglvl, equation->b, equation->x, &error);

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];
    for (int iDim = 0; iDim < 3; iDim++) {
      if (!node->bc[iDim]) {
        node->u[iDim] = equation->x[3 * iNode + iDim];
      } else {
        node->u[iDim] = scale_u * node->bc_u[iDim];
      }
    }
  }
}
#else
#include "dmumps_c.h"

#define ICNTL(I) icntl[(I) - 1]

static void SetDirichletBCForMatrix(FEMEQUATION_t *equation) {
  int n = equation->n;
  double *a = equation->a;
  double *bc = equation->bc;

  for (int i = 0; i < n; i++) {
    int iDofID = equation->row[i];
    int jDofID = equation->col[i];
    if (bc[iDofID] < 0.5) {
      if (iDofID == jDofID) {
        a[i] = 1.0;
      } else {
        a[i] = 0.0;
      }
    } else if (bc[jDofID] < 0.5) {
      a[i] = 0.0;
    }
  }
}

void CheckErrorInMUMPS(char *phase, int info) {
  if (info != 0) {
    fprintf (stdout, "[MSG] Error in MUMPS happens at %s phase.\n", phase);
    exit(1);
  }
}

void FEM_Factorization(FEM_t *fem) {
  int nNodes = fem->nNodes;
  FEMEQUATION_t *equation = &fem->equation;
  DMUMPS_STRUC_C *id = &equation->id;

  if (fem->nElements == 0) return;

  fprintf(stdout, "[LOG] Compute the inverse of the stiffness matrix\n");
  fflush(stdout);

  SetDirichletBCForMatrix(equation);

  for (int i = 0; i < equation->n; i++) {
    equation->row[i] += 1;
    equation->col[i] += 1;
  }

  id->job = -1;
  id->par = 1;
  id->sym = 2;
  id->comm_fortran = -987654;
  dmumps_c(id);

  CheckErrorInMUMPS("Initialization", id->INFO(1));
  id->ICNTL(1) = -1;
  id->ICNTL(2) = -1;
  id->ICNTL(3) = -1;
  id->ICNTL(4) = 0;

  id->job = 4;
  id->n = 3 * nNodes;
  id->nnz = equation->n;
  id->irn = equation->row;
  id->jcn = equation->col;
  id->a = equation->a;
  id->rhs = equation->b;
  dmumps_c(id);
  CheckErrorInMUMPS("Factorization", id->INFO(1));

  for (int i = 0; i < equation->n; i++) {
    equation->row[i] -= 1;
    equation->col[i] -= 1;
  }
}

static void SetDirichletBCForRHS(FEMEQUATION_t *equation) {
  int n = equation->n;
  double *a = equation->a0;
  double *x = equation->x;
  double *b = equation->b;
  double *bc = equation->bc;

  for (int i = 0; i < n; i++) {
    int iDofID = equation->row[i];
    int jDofID = equation->col[i];

    if (bc[iDofID] < 0.5) {
      if (iDofID != jDofID) {
        b[jDofID] -= a[i] * x[iDofID];
      }
    } else if (bc[jDofID] < 0.5) {
      b[iDofID] -= a[i] * x[jDofID];
    }
  }
}

static void MakeUnknown(FEM_t *fem) {
  double *x = fem->equation.x;
  double scale_u = fem->equation.scale_u;

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    for (int iDim = 0; iDim < 3; iDim++) {
      x[3 * iNode + iDim] = scale_u * node->bc_u[iDim];
    }
  }
}

static void MakeRHS(FEM_t *fem) {
  double *b = fem->equation.b;
  double scale_f = fem->equation.scale_f;

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    for (int iDim = 0; iDim < 3; iDim++) {
      b[3 * iNode + iDim] = node->f[iDim] + scale_f * node->bc_f[iDim];
    }
  }

  SetDirichletBCForRHS(&fem->equation);
}

void FEM_ComputeSolution(FEM_t *fem) {
  FEMEQUATION_t *equation = &fem->equation;
  DMUMPS_STRUC_C *id = &equation->id;
  double scale_u = equation->scale_u;

  if (fem->nElements == 0) return;

  MakeUnknown(fem);
  MakeRHS(fem);

  for (int i = 0; i < equation->n; i++) {
    equation->row[i] += 1;
    equation->col[i] += 1;
  }

  id->job = 3;
  dmumps_c(id);
  CheckErrorInMUMPS("Solution", id->INFO(1));

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];
    for (int iDim = 0; iDim < 3; iDim++) {
      if (!node->bc[iDim]) {
        node->u[iDim] = equation->b[3 * iNode + iDim];
      } else {
        node->u[iDim] = scale_u * node->bc_u[iDim];
      }
    }
  }
  for (int i = 0; i < equation->n; i++) {
    equation->row[i] -= 1;
    equation->col[i] -= 1;
  }
}
#endif
#else
void FEM_Factorization(FEM_t *fem) { return; }

void FEM_ComputeSolution(FEM_t *fem) { return; }
#endif
