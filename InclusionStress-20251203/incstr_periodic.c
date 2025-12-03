#include <stdio.h>
#include <stdlib.h>

#include "incstr_struct.h"

void INCSTR_MakePeriodicImages(MATERIAL_t *material) {
  material->offset = (double **)malloc(27 * sizeof(double *));
  for (int i = 0; i < 27; i++) {
    material->offset[i] = (double *)malloc(3 * sizeof(double));
  }

  material->nImages = 1;
  for (int iDim = 0; iDim < 3; iDim++) {
    material->offset[0][iDim] = 0.0;
  }

  if (material->periodic[0]) {
    material->nImages = 3;
    material->offset[1][0] = material->size[0];
    material->offset[1][1] = 0.0;
    material->offset[1][2] = 0.0;
    material->offset[2][0] = -1.0 * material->size[0];
    material->offset[2][1] = 0.0;
    material->offset[2][2] = 0.0;

    if (material->periodic[1]) {
      material->nImages = 9;
      material->offset[3][0] = 0.0;
      material->offset[3][1] = material->size[1];
      material->offset[3][2] = 0.0;
      material->offset[4][0] = 0.0;
      material->offset[4][1] = -1.0 * material->size[1];
      material->offset[4][2] = 0.0;

      material->offset[5][0] = material->size[0];
      material->offset[5][1] = material->size[1];
      material->offset[5][2] = 0.0;
      material->offset[6][0] = -1.0 * material->size[0];
      material->offset[6][1] = material->size[1];
      material->offset[6][2] = 0.0;
      material->offset[7][0] = material->size[0];
      material->offset[7][1] = -1.0 * material->size[1];
      material->offset[7][2] = 0.0;
      material->offset[8][0] = -1.0 * material->size[0];
      material->offset[8][1] = -1.0 * material->size[1];
      material->offset[8][2] = 0.0;

      if (material->periodic[2]) {
        material->nImages = 27;
        material->offset[9][0] = 0.0;
        material->offset[9][1] = 0.0;
        material->offset[9][2] = material->size[2];
        material->offset[10][0] = 0.0;
        material->offset[10][1] = 0.0;
        material->offset[10][2] = -1.0 * material->size[2];

        material->offset[11][0] = material->size[0];
        material->offset[11][1] = 0.0;
        material->offset[11][2] = material->size[2];
        material->offset[12][0] = material->size[0];
        material->offset[12][1] = 0.0;
        material->offset[12][2] = -1.0 * material->size[2];
        material->offset[13][0] = -1.0 * material->size[0];
        material->offset[13][1] = 0.0;
        material->offset[13][2] = material->size[2];
        material->offset[14][0] = -1.0 * material->size[0];
        material->offset[14][1] = 0.0;
        material->offset[14][2] = -1.0 * material->size[2];

        material->offset[15][0] = 0.0;
        material->offset[15][1] = material->size[1];
        material->offset[15][2] = material->size[2];
        material->offset[16][0] = 0.0;
        material->offset[16][1] = material->size[1];
        material->offset[16][2] = -1.0 * material->size[2];
        material->offset[17][0] = 0.0;
        material->offset[17][1] = -1.0 * material->size[1];
        material->offset[17][2] = material->size[2];
        material->offset[18][0] = 0.0;
        material->offset[18][1] = -1.0 * material->size[1];
        material->offset[18][2] = -1.0 * material->size[2];

        material->offset[19][0] = material->size[0];
        material->offset[19][1] = material->size[1];
        material->offset[19][2] = material->size[2];
        material->offset[20][0] = -1.0 * material->size[0];
        material->offset[20][1] = material->size[1];
        material->offset[20][2] = material->size[2];
        material->offset[21][0] = -1.0 * material->size[0];
        material->offset[21][1] = -1.0 * material->size[1];
        material->offset[21][2] = material->size[2];
        material->offset[22][0] = material->size[0];
        material->offset[22][1] = -1.0 * material->size[1];
        material->offset[22][2] = material->size[2];
        material->offset[23][0] = material->size[0];
        material->offset[23][1] = material->size[1];
        material->offset[23][2] = -1.0 * material->size[2];
        material->offset[24][0] = -1.0 * material->size[0];
        material->offset[24][1] = material->size[1];
        material->offset[24][2] = -1.0 * material->size[2];
        material->offset[25][0] = -1.0 * material->size[0];
        material->offset[25][1] = -1.0 * material->size[1];
        material->offset[25][2] = -1.0 * material->size[2];
        material->offset[26][0] = material->size[0];
        material->offset[26][1] = -1.0 * material->size[1];
        material->offset[26][2] = -1.0 * material->size[2];
      }
    } else {
      if (material->periodic[2]) {
        material->nImages = 9;
        material->offset[3][0] = 0.0;
        material->offset[3][1] = 0.0;
        material->offset[3][2] = material->size[2];
        material->offset[4][0] = 0.0;
        material->offset[4][1] = 0.0;
        material->offset[4][2] = -1.0 * material->size[2];

        material->offset[5][0] = material->size[0];
        material->offset[5][1] = 0.0;
        material->offset[5][2] = material->size[2];
        material->offset[6][0] = material->size[0];
        material->offset[6][1] = 0.0;
        material->offset[6][2] = -1.0 * material->size[2];
        material->offset[7][0] = -1.0 * material->size[0];
        material->offset[7][1] = 0.0;
        material->offset[7][2] = material->size[2];
        material->offset[8][0] = -1.0 * material->size[0];
        material->offset[8][1] = 0.0;
        material->offset[8][2] = -1.0 * material->size[2];
      }
    }
  } else {
    if (material->periodic[1]) {
      material->nImages = 3;
      material->offset[1][0] = 0.0;
      material->offset[1][1] = material->size[1];
      material->offset[1][2] = 0.0;
      material->offset[2][0] = 0.0;
      material->offset[2][1] = -1.0 * material->size[1];
      material->offset[2][2] = 0.0;

      if (material->periodic[2]) {
        material->nImages = 9;
        material->offset[3][0] = 0.0;
        material->offset[3][1] = 0.0;
        material->offset[3][2] = material->size[2];
        material->offset[4][0] = 0.0;
        material->offset[4][1] = 0.0;
        material->offset[4][2] = -1.0 * material->size[2];

        material->offset[5][0] = 0.0;
        material->offset[5][1] = material->size[1];
        material->offset[5][2] = material->size[2];
        material->offset[6][0] = 0.0;
        material->offset[6][1] = material->size[1];
        material->offset[6][2] = -1.0 * material->size[2];
        material->offset[7][0] = 0.0;
        material->offset[7][1] = -1.0 * material->size[1];
        material->offset[7][2] = material->size[2];
        material->offset[8][0] = 0.0;
        material->offset[8][1] = -1.0 * material->size[1];
        material->offset[8][2] = -1.0 * material->size[2];
      }
    } else {
      if (material->periodic[2]) {
        material->nImages = 3;
        material->offset[1][0] = 0.0;
        material->offset[1][1] = 0.0;
        material->offset[1][2] = material->size[2];
        material->offset[2][0] = 0.0;
        material->offset[2][1] = 0.0;
        material->offset[2][2] = -1.0 * material->size[2];
      }
    }
  }
}