#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include "dd_defs.h"
#include "dd_struct.h"

static FILE *OpenPlot(DD_t *dd) {
  FILE *gfp = popen(GNUPLOT_PATH, "w");
  int font_size = 15;
  double point_size = 0.5;

  fprintf(gfp, "set border 4095\n");
  fprintf(gfp, "set xrange [%e:%e]\n", 0.0, dd->bc.size[0]);
  fprintf(gfp, "set yrange [%e:%e]\n", 0.0, dd->bc.size[1]);
  fprintf(gfp, "set zrange [%e:%e]\n", 0.0, dd->bc.size[2]);
  fprintf(gfp, "set xlabel font \"Arial, %d\"\n", font_size);
  fprintf(gfp, "set ylabel font \"Arial, %d\"\n", font_size);
  fprintf(gfp, "set zlabel font \"Arial, %d\"\n", font_size);
  fprintf(gfp, "set xlabel \"X\"\n");
  fprintf(gfp, "set ylabel \"Y\"\n");
  fprintf(gfp, "set zlabel \"Z\"\n");
  /*
  fprintf(gfp, "set noxtics\n");
  fprintf(gfp, "set noytics\n");
  fprintf(gfp, "set noztics\n");
  */
  fprintf(gfp, "set ticslevel 0\n");
  fprintf(gfp, "set pointsize %f\n", point_size);

  return gfp;
}

static int CheckOutputElements(DD_t *dd) {
  int nElements = 0;
  double *size = dd->bc.size;
  int *periodic = dd->bc.periodic;

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];
    int nodeID0 = element->nodeID[0];
    int nodeID1 = element->nodeID[1];
    NODE_t *node0 = &dd->node[nodeID0];
    NODE_t *node1 = &dd->node[nodeID1];

    // initial value: output
    element->iTmp = 1;

    for (int iDim = 0; iDim < 3; iDim++) {
      double dx = node1->x[iDim] - node0->x[iDim];

      if (periodic[iDim] && fabs(dx) > 0.5 * size[iDim]) {
        element->iTmp = 0;
        break;
      }
    }

    if (element->iTmp) {
      nElements += 1;
    }
  }

  return nElements;
}

void DD_PlotDislocations(DD_t *dd) {
  FILE *fp = fopen("elem.plt", "w");
  static int first_time = 1;
  static FILE *gfp;

  if (first_time) {
    gfp = OpenPlot(dd);
    first_time = 0;
  }

  // Need some time to process
  usleep(100000);

  CheckOutputElements(dd);
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (!element->iTmp) continue;

    for (int iNode = 0; iNode < 2; iNode++) {
      int nodeID = element->nodeID[iNode];
      NODE_t *node = &dd->node[nodeID];

      for (int iDim = 0; iDim < 3; iDim++) {
        fprintf(fp, "%e ", node->x[iDim]);
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fprintf(fp, "\n");
  }
  fclose(fp);

  fprintf(gfp, "splot \"./elem.plt\" notitle w lp lt 7 lc 7\n");
  fflush(gfp);
}
