#include <stdio.h>

#include "dd_element_handle.h"
#include "dd_struct.h"
#include "fem_struct.h"

static double SimulationVolume(double *s) {
  double volume = 1.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    volume *= s[iDim];
  }
  return volume;
}

void DD_FEM_WriteDislocationDensity(DD_t *dd, FEM_t *fem, char *directory) {
  char file_name[256];
  FILE *fp;
  double volume = 0.0;
  double total_length = 0.0;
  double density = 0.0;

  sprintf(file_name, "%s/density.out", directory);
  if (dd->output.id == 1) {
    fp = fopen(file_name, "w");
  } else {
    fp = fopen(file_name, "a");
  }

  if (fem->nElements == 0) {
    volume = SimulationVolume(dd->bc.size);
  } else {
    volume = fem->volume;
  }

  DD_TangentVectors(dd);
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    total_length += dd->element[iElement].length;
  }
  density = total_length / volume;
  fprintf(fp, "%e %e\n", dd->step.t, density);

  fclose(fp);
}
