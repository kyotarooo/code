#include <stdio.h>

#include "fem_struct.h"

void FEM_WriteNodes(FEM_t *fem, char *directory) {
  char file_name[256];
  FILE *fp;

  sprintf(file_name, "%s/node.out", directory);
  fp = fopen(file_name, "w");

  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    fprintf(fp, "%d", iNode + 1);

    for (int iDim = 0; iDim < 3; iDim++) {
      fprintf(fp, " %e", node->u[iDim]);
    }
    for (int iDim = 0; iDim < 6; iDim++) {
      fprintf(fp, " %e", node->stress[iDim]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

void FEM_WriteMesh(FEM_t *fem, char *directory) {
  char file_name[256];
  FILE *fp;
  int node_index[20] = {0, 1, 2,  3,  12, 13, 14, 15, 4,  5,
                        6, 7, 16, 17, 18, 19, 8,  9,  10, 11};

  if (fem->nElements == 0) return;

  sprintf(file_name, "%s/mesh.vtk", directory);
  fp = fopen(file_name, "w");

  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "FEM Mesh\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fp, "POINTS %d float\n", fem->nNodes);
  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    for (int iDim = 0; iDim < 3; iDim++) {
      fprintf(fp, "%e ", node->x[iDim]);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "CELLS %d %d\n", fem->nElements, fem->nElements * 21);
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];

    fprintf(fp, "20 ");
    for (int iNode = 0; iNode < 20; iNode++) {
      fprintf(fp, "%d ", element->nodeID[node_index[iNode]]);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "CELL_TYPES %d\n", fem->nElements);
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    fprintf(fp, "25\n");
  }

  fclose(fp);
}

void FEM_WriteStress(int stressID, FEM_t *fem, char *directory) {
  char file_name[256];
  FILE *fp;
  int node_index[20] = {0, 1, 2,  3,  12, 13, 14, 15, 4,  5,
                        6, 7, 16, 17, 18, 19, 8,  9,  10, 11};

  if (fem->nElements == 0) return;

  sprintf(file_name, "%s/stress.vtk", directory);
  fp = fopen(file_name, "w");

  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "FEM Mesh\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fp, "POINTS %d float\n", fem->nNodes);
  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    FEMNODE_t *node = &fem->node[iNode];

    for (int iDim = 0; iDim < 3; iDim++) {
      fprintf(fp, "%e ", node->x[iDim]);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "CELLS %d %d\n", fem->nElements, fem->nElements * 21);
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];

    fprintf(fp, "20 ");
    for (int iNode = 0; iNode < 20; iNode++) {
      fprintf(fp, "%d ", element->nodeID[node_index[iNode]]);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "CELL_TYPES %d\n", fem->nElements);
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    fprintf(fp, "25\n");
  }

  fprintf (fp, "POINT_DATA %d\n", fem->nNodes);
  fprintf (fp, "SCALARS stress float\n");
  fprintf (fp, "LOOKUP_TABLE default\n");
  for(int iNode = 0; iNode < fem->nNodes; iNode++) {
    fprintf (fp, "%e\n", fem->node[iNode].stress[stressID]);
  }

  fclose(fp);
}