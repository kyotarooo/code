#include <stdio.h>

#include "dd_dislocation_stress.h"
#include "dd_grid_stress.h"
#include "dd_pk_force.h"
#include "dd_struct.h"
#include "fem_stress.h"
#include "fem_struct.h"

void DD_FEM_ElasticInteractionForce(double x[3], double b[3], double t[3],
                                    DD_t *dd, FEM_t *fem, double f[3]) {
  double applied_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double dislocation_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double grid_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double correction_stress[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double stress[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  // Dislocation-dislocation interaction
  DD_DislocationStress(x, dd, dislocation_stress);

  // Dislocation-initial stress interaction
  DD_GridStress(x, dd, grid_stress);

  // Either applied stress or correction stress
  if (fem->nElements == 0) {
    // If we have FE mesh, applied stress is ignored (NOT APPLIED).
    for (int iDim = 0; iDim < 3; iDim++) {
      for (int jDim = 0; jDim < 3; jDim++) {
        applied_stress[iDim][jDim] = dd->bc.stress[iDim][jDim];
      }
    }
  } else {
    // If we do not have FE mesh, applied stress is APPLIED.
    // Correction stress is calculated to consider the dislocation-surface
    // interaction
    FEM_CorrectionStress(x, fem, correction_stress);
  }

  // Sum all stress including applied stress
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      stress[iDim][jDim] +=
          applied_stress[iDim][jDim] + dislocation_stress[iDim][jDim] +
          grid_stress[iDim][jDim] + correction_stress[iDim][jDim];
    }
  }

  // Convert stress to force on dislocation
  DD_PeachKoehlerForce(stress, b, t, f);
}
