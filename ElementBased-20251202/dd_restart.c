#include <stdio.h>
#include <stdlib.h>

#include "dd_defs.h"
#include "dd_struct.h"

static void ReadIntegers(int *v, int n, FILE *fp) {
  fread(v, sizeof(int), n, fp);
}

static void ReadReals(double *v, int n, FILE *fp) {
  fread(v, sizeof(double), n, fp);
}

static void WriteIntegers(int *v, int n, FILE *fp) {
  fwrite(v, sizeof(int), n, fp);
}

static void WriteReals(double *v, int n, FILE *fp) {
  fwrite(v, sizeof(double), n, fp);
}

void DD_Restart(int *stepID, DD_t *dd, char *directory, int mode) {
  char file_name[256];
  FILE *fp;
  void (*integers)(int *v, int n, FILE *fp);
  void (*reals)(double *v, int n, FILE *fp);

  sprintf(file_name, "%s/restart.dat", directory);
  switch (mode) {
    case READ_RESTART_DATA:
      fp = fopen(file_name, "r");
      integers = ReadIntegers;
      reals = ReadReals;
      break;
    case WRITE_RESTART_DATA:
      fp = fopen(file_name, "w");
      integers = WriteIntegers;
      reals = WriteReals;
      break;
  }

  // if no restart data
  if (fp == NULL) {
    return;
  }

  integers(stepID, 1, fp);

  integers(&dd->nNodes, 1, fp);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];

    integers(&node->type, 1, fp);
    reals(node->x, 3, fp);
    reals(node->x0, 3, fp);
    reals(node->v, 3, fp);
    reals(node->direction, 3, fp);
  }

  integers(&dd->nElements, 1, fp);
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    integers(element->nodeID, 2, fp);
    reals(element->burgers, 3, fp);
    reals(element->slip, 3, fp);
    reals(&element->length, 1, fp);
    reals(element->tangent, 3, fp);

    for (int iIntegral = 0; iIntegral < N_DD_INTEGRALS; iIntegral++) {
      INTEGRAL_t *integral = &element->integral[iIntegral];

      reals(integral->core_x, dd->core.nDislocations, fp);
    }
  }

  reals(&dd->step.t, 1, fp);
  reals(&dd->step.dt, 1, fp);

  for (int iDim = 0; iDim < 3; iDim++) {
    reals(dd->bc.stress[iDim], 3, fp);
  }

  for (int iDim = 0; iDim < 3; iDim++) {
    reals(dd->mechanical_behavior.strain[iDim], 3, fp);
    reals(dd->mechanical_behavior.plastic_strain[iDim], 3, fp);
  }

  integers(&dd->output.id, 1, fp);

  fclose(fp);
}
