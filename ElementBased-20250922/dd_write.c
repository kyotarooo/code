#include <math.h>
#include <stdint.h>
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

void DD_WriteTime(DD_t *dd, char *directory) {
  char file_name[256];
  FILE *fp;

  // Time data file
  sprintf(file_name, "%s/time.out", directory);
  fp = fopen(file_name, "a");
  fprintf(fp, "%d %e\n", dd->output.id, dd->step.t);
  fclose(fp);
}

#ifdef _OUTPUT_VTU_
void DD_WriteDislocations(DD_t *dd, char *directory) {
  int nElements = CheckOutputElements(dd);
  char file_name[256];
  int id;
  double t;
  int offset = 0;
  FILE *fp, *tfp;

  // PVD file
  sprintf(file_name, "%s/time.out", directory);
  tfp = fopen(file_name, "r");
  sprintf(file_name, "%s/elem.pvd", directory);
  fp = fopen(file_name, "w");
  fprintf(fp, "<VTKFile type=\"Collection\">\n");
  fprintf(fp, "  <Collection>\n");
  while (fscanf(tfp, "%d%lf", &id, &t) != EOF) {
    fprintf(
        fp,
        "    <DataSet part=\"0\" timestep=\"%e\" file=\"elem_%04d.vtu\" />\n",
        t, id);
  }
  fprintf(fp, "  </Collection>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(tfp);
  fclose(fp);

  // VTU file
  sprintf(file_name, "%s/elem_%04d.vtu", directory, dd->output.id);
  fp = fopen(file_name, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          dd->nNodes, nElements);
  fprintf(fp, "      <Points>\n");
  fprintf(fp,
          "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" "
          "format=\"ascii\">\n");
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    double *x = dd->node[iNode].x;
    fprintf(fp, "          ");
    for (int iDim = 0; iDim < 3; iDim++) {
      fprintf(fp, "%e ", x[iDim]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");
  fprintf(fp, "      <Cells>\n");
  fprintf(fp,
          "        <DataArray type=\"Int32\" Name=\"connectivity\" "
          "format=\"ascii\">\n");
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (!element->iTmp) continue;

    fprintf(fp, "          ");
    for (int iNode = 0; iNode < 2; iNode++) {
      fprintf(fp, "%d ", element->nodeID[iNode]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(
      fp,
      "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (!element->iTmp) continue;

    offset += 2;
    fprintf(fp, "          %d\n", offset);
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(
      fp,
      "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (!element->iTmp) continue;

    fprintf(fp, "          3\n");
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp,
          "        <DataArray type=\"Int32\" Name=\"NodeType\" "
          "format=\"ascii\">\n");
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    fprintf(fp, "          %d\n", dd->node[iNode].type);
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </PointData>\n");
  fprintf(fp, "      <CellData Scalars=\"scalars\">\n");
  fprintf(fp,
          "        <DataArray type=\"Float32\" Name=\"BurgersVector\" "
          "NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (!element->iTmp) continue;

    fprintf(fp, "          ");
    for (int iDim = 0; iDim < 3; iDim++) {
      fprintf(fp, "%e ", element->burgers[iDim]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp,
          "        <DataArray type=\"Float32\" Name=\"SlipPlaneNormalVector\" "
          "NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (!element->iTmp) continue;

    fprintf(fp, "          ");
    for (int iDim = 0; iDim < 3; iDim++) {
      fprintf(fp, "%e ", element->slip[iDim]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </CellData>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");

  fclose(fp);
}

#else

void DD_WriteDislocations(DD_t *dd, char *directory) {
  int nElements = CheckOutputElements(dd);
  char file_name[256];
  FILE *fp;

  sprintf(file_name, "%s/elem_%04d.vtk", directory, dd->output.id);
  fp = fopen(file_name, "w");

  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "%e\n", dd->step.t);
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fp, "POINTS %d float\n", dd->nNodes);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    double *x = dd->node[iNode].x;

    for (int iDim = 0; iDim < 3; iDim++) {
      fprintf(fp, "%e ", x[iDim]);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "CELLS %d %d\n", nElements, nElements * 3);
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (!element->iTmp) continue;

    fprintf(fp, "2 ");
    for (int iNode = 0; iNode < 2; iNode++) {
      fprintf(fp, "%d ", element->nodeID[iNode]);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "CELL_TYPES %d\n", nElements);
  for (int iElement = 0; iElement < nElements; iElement++) {
    fprintf(fp, "3\n");
  }

  fprintf(fp, "POINT_DATA %d\n", dd->nNodes);
  fprintf(fp, "SCALARS NodeType float\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    fprintf(fp, "%d\n", dd->node[iNode].type);
  }
  fprintf(fp, "CELL_DATA %d\n", nElements);
  fprintf(fp, "VECTORS BurgersVector float\n");
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (!element->iTmp) continue;

    for (int iDim = 0; iDim < 3; iDim++) {
      fprintf(fp, "%e ", element->burgers[iDim]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "VECTORS SlipPlaneNormalVector float\n");
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    if (!element->iTmp) continue;

    for (int iDim = 0; iDim < 3; iDim++) {
      fprintf(fp, "%e ", element->slip[iDim]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}
#endif

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
