#include <math.h>
#include <stdio.h>
#include <time.h>

#include "dd_struct.h"

void DD_WriteLog(int stepID, DD_t *dd) {
  int i = dd->mechanical_behavior.index[0];
  int j = dd->mechanical_behavior.index[1];
  time_t current_time = time(NULL);
  char *current_time_string = ctime(&current_time);

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "- Current time: %s", current_time_string);
  fprintf(stdout, "----------------------------------------\n");
  fprintf(stdout, "- Step: %d Time: %e\n", stepID, dd->step.t);
  fprintf(stdout, "- Nodes: %d Elements: %d\n", dd->nNodes, dd->nElements);
  fprintf(stdout, "- Stress:\n");
  for (int iDim = 0; iDim < 3; iDim++) {
    fprintf(stdout, "  ");
    for (int jDim = 0; jDim < 3; jDim++) {
      fprintf(stdout, "%e ", dd->bc.stress[iDim][jDim]);
    }
    fprintf(stdout, "\n");
  }

  if (i != -1 && j != -1) {
    fprintf(stdout, "- Total strain [%d %d]: %e\n", i, j,
            dd->mechanical_behavior.strain[i][j]);
    fprintf(stdout, "- Plastic strain [%d %d]: %e\n", i, j,
            dd->mechanical_behavior.plastic_strain[i][j]);
  }

  fflush(stdout);
}

static int DeleteElementID(DD_t *dd) {
  double *size = dd->bc.size;
  int *periodic = dd->bc.periodic;

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];
    int nodeID0 = element->nodeID[0];
    int nodeID1 = element->nodeID[1];
    NODE_t *node0 = &dd->node[nodeID0];
    NODE_t *node1 = &dd->node[nodeID1];

    for (int iDim = 0; iDim < 3; iDim++) {
      double dx = node1->x[iDim] - node0->x[iDim];

      if (periodic[iDim] && fabs(dx) > 0.5 * size[iDim]) {
        return iElement;
      }
    }
  }

  return -1;
}

typedef enum { Glide, Line, Plane } Direction;

void DD_WriteDislocations(DD_t *dd, char *directory) {
  Direction glide_direction = Glide;
  Direction line_direction = Line;
  Direction plane_direction = Plane;
  int delete_elementID = DeleteElementID(dd);
  ELEMENT_t *delete_element = &dd->element[delete_elementID];
  int nodeID0 = delete_element->nodeID[0];
  int nodeID1 = delete_element->nodeID[1];
  NODE_t *node0 = &dd->node[nodeID0];
  NODE_t *node1 = &dd->node[nodeID1];
  double g = 0.5 * (node0->x[pbc_dir] + node1->x[pbc_dir]);
  char file_name[256];
  FILE *fp;

  sprintf(file_name, "%s/elem_%04d.vtk", directory, dd->output.id);
  fp = fopen(file_name, "w");

  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "%e\n", dd->step.t);
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fp, "POINTS %d float\n", dd->nNodes + 2);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];

    for (int iDim = 0; iDim < 3; iDim++) {
      if (iDim == glide_direction) {
        fprintf(fp, "%e ", g);
      } else if (iDim == line_direction) {
        fprintf(fp, "%e ", 0.0);
      } else {
        fprintf(fp, "%e ", node0->x[iDim]);
      }
    }
    fprintf(fp, "\n");
    for (int iDim = 0; iDim < 3; iDim++) {
      if (iDim == glide_direction) {
        fprintf(fp, "%e ", g);
      } else if (iDim == line_direction) {
        fprintf(fp, "%e ", size[iDim]);
      } else {
        fprintf(fp, "%e ", node0->x[iDim]);
      }
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "CELLS %d %d\n", (dd->nElements - 1 + 2),
          (dd->nElements - 1 + 2) * 3);
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (iElement == delete_elementID) {
      if (node0->x[line_direction] < 0.5 * size[line_direction]) {
        fprintf(fp, "2 %d %d\n", nodeID0, dd->nNodes);
        fprintf(fp, "2 %d %d\n", dd->nNodes + 1, nodeID1);
      } else {
        fprintf(fp, "2 %d %d\n", nodeID0, dd->nNodes + 1);
        fprintf(fp, "2 %d %d\n", dd->nNodes, nodeID1);
      }
    } else {
      fprintf(fp, "2 ");
      for (int iNode = 0; iNode < 2; iNode++) {
        fprintf(fp, "%d ", element->nodeID[iNode]);
      }
      fprintf(fp, "\n");
    }
  }

  fprintf(fp, "CELL_TYPES %d\n", dd->nElements - 1 + 2);
  for (int iElement = 0; iElement < dd->nElements - 1 + 2; iElement++) {
    fprintf(fp, "3\n");
  }

  /*
    fprintf(fp, "CELL_DATA %d\n", nElements);
    fprintf(fp, "VECTORS Burgers float\n");
    for (int iElement = 0; iElement < dd->nElements; iElement++) {
      ELEMENT_t *element = &dd->element[iElement];

      if (!element->iTmp) continue;

      for (int iDim = 0; iDim < 3; iDim++) {
        fprintf(fp, "%e ", element->burgers[iDim]);
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "VECTORS SlipPlane float\n");
    for (int iElement = 0; iElement < dd->nElements; iElement++) {
      ELEMENT_t *element = &dd->element[iElement];

      if (!element->iTmp) continue;

      for (int iDim = 0; iDim < 3; iDim++) {
        fprintf(fp, "%e ", element->slip[iDim]);
      }
      fprintf(fp, "\n");
    }
  */
  fclose(fp);
}

void DD_WriteMechanicalBehavior(DD_t *dd, char *directory) {
  char file_name[256];
  FILE *fp;
  int i = dd->mechanical_behavior.index[0];
  int j = dd->mechanical_behavior.index[1];

  if (i == -1 || j == -1) return;

  sprintf(file_name, "%s/mechanical_behavior.out", directory);
  if (dd->output.id == 1) {
    fp = fopen(file_name, "w");
  } else {
    fp = fopen(file_name, "a");
  }

  fprintf(fp, "%e %e %e %e\n", dd->step.t, dd->mechanical_behavior.strain[i][j],
          dd->mechanical_behavior.plastic_strain[i][j], dd->bc.stress[i][j]);

  fclose(fp);
}

void DD_WriteSimulationVolume(DD_t *dd, char *directory) {
  char file_name[256];
  FILE *fp;
  double x0 = 0.0;
  double x1 = dd->bc.size[0];
  double y0 = 0.0;
  double y1 = dd->bc.size[1];
  double z0 = 0.0;
  double z1 = dd->bc.size[2];

  sprintf(file_name, "%s/volume.vtk", directory);
  fp = fopen(file_name, "w");

  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "%e\n", dd->step.t);
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fp, "POINTS %d float\n", 8);
  fprintf(fp, "%e %e %e\n", x0, y0, z0);
  fprintf(fp, "%e %e %e\n", x1, y0, z0);
  fprintf(fp, "%e %e %e\n", x1, y1, z0);
  fprintf(fp, "%e %e %e\n", x0, y1, z0);
  fprintf(fp, "%e %e %e\n", x0, y0, z1);
  fprintf(fp, "%e %e %e\n", x1, y0, z1);
  fprintf(fp, "%e %e %e\n", x1, y1, z1);
  fprintf(fp, "%e %e %e\n", x0, y1, z1);

  fprintf(fp, "CELLS %d %d\n", 1, 9);
  fprintf(fp, "8 0 1 2 3 4 5 6 7\n");

  fprintf(fp, "CELL_TYPES 1\n");
  fprintf(fp, "12\n");

  fclose(fp);
}
