#include <math.h>
#include <stdlib.h>

#include "fem_struct.h"

void FEM_FindSurface(FEM_t *fem) {
  int **node_to_element;
  int surface_node_index[6][8] = {
      {3, 2, 1, 0, 6, 5, 4, 7},     {0, 1, 13, 12, 4, 9, 16, 8},
      {1, 2, 14, 13, 5, 10, 17, 9}, {2, 3, 15, 14, 6, 11, 18, 10},
      {3, 0, 12, 15, 7, 11, 19, 8}, {12, 13, 14, 15, 16, 17, 18, 19}};

  // Allocate memory
  node_to_element = (int **)malloc(fem->nNodes * sizeof(int *));

  // Initialize
  for (int iNode = 0; iNode < fem->nNodes; iNode++) fem->node[iNode].iTmp = 0;

  // Count the number of elements using a node
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];
    for (int iNode = 0; iNode < 20; iNode++) {
      int nodeID = element->nodeID[iNode];

      fem->node[nodeID].iTmp++;
    }
  }

  // Allocate memory
  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    node_to_element[iNode] = (int *)malloc(fem->node[iNode].iTmp * sizeof(int));
  }

  // Initialize
  for (int iNode = 0; iNode < fem->nNodes; iNode++) fem->node[iNode].iTmp = 0;

  // Count the number of elements using a node
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];
    for (int iNode = 0; iNode < 20; iNode++) {
      int nodeID = element->nodeID[iNode];

      node_to_element[nodeID][fem->node[nodeID].iTmp] = iElement;
      fem->node[nodeID].iTmp++;
    }
  }

  // Search the surface
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];

    // Initialize
    element->nSurfaces = 6;

    for (int iSurface = 0; iSurface < 6; iSurface++) {
      int nodeID0 = element->nodeID[surface_node_index[iSurface][0]];
      int nodeID1 = element->nodeID[surface_node_index[iSurface][2]];
      int nElements0 = fem->node[nodeID0].iTmp;
      int nElements1 = fem->node[nodeID1].iTmp;

      // Initialize
      element->surface[iSurface] = 1;

      for (int jElement = 0; jElement < nElements0; jElement++) {
        int elementID0 = node_to_element[nodeID0][jElement];

        if (elementID0 == iElement) continue;

        for (int kElement = 0; kElement < nElements1; kElement++) {
          int elementID1 = node_to_element[nodeID1][kElement];
          // The surface is commonly used -> it is not surface
          if (elementID0 == elementID1) {
            element->surface[iSurface] = 0;
            element->nSurfaces -= 1;
            break;
          }
        }

        if (!element->surface[iSurface]) break;
      }
    }
  }

  // Release memory
  for (int iNode = 0; iNode < fem->nNodes; iNode++) {
    free(node_to_element[iNode]);
  }
  free(node_to_element);

  // Center of element, center of element face and its unit normal vector
  for (int iElement = 0; iElement < fem->nElements; iElement++) {
    FEMELEMENT_t *element = &fem->element[iElement];

    element->surface_nx = (double **)malloc(6 * sizeof(double *));
    element->surface_x = (double **)malloc(6 * sizeof(double *));
    element->x = (double *)malloc(3 * sizeof(double));

    for (int iDim = 0; iDim < 3; iDim++) {
      double x = 0.0;
      for (int iNode = 0; iNode < 20; iNode++) {
        int nodeID = element->nodeID[iNode];

        x += fem->node[nodeID].x[iDim];
      }

      element->x[iDim] = x / 20.0;
    }

    for (int iSurface = 0; iSurface < 6; iSurface++) {
      element->surface_nx[iSurface] = (double *)malloc(3 * sizeof(double));
      element->surface_x[iSurface] = (double *)malloc(3 * sizeof(double));

      for (int iDim = 0; iDim < 3; iDim++) {
        double x = 0.0;
        for (int iNode = 0; iNode < 8; iNode++) {
          int nodeID = element->nodeID[surface_node_index[iSurface][iNode]];
          FEMNODE_t *node = &fem->node[nodeID];

          x += node->x[iDim];
        }

        element->surface_x[iSurface][iDim] = x / 8.0;
      }

      {
        int nodeID0 = element->nodeID[surface_node_index[iSurface][0]];
        int nodeID1 = element->nodeID[surface_node_index[iSurface][1]];
        int nodeID2 = element->nodeID[surface_node_index[iSurface][2]];
        FEMNODE_t *node0 = &fem->node[nodeID0];
        FEMNODE_t *node1 = &fem->node[nodeID1];
        FEMNODE_t *node2 = &fem->node[nodeID2];
        double rx0 = node1->x[0] - node0->x[0];
        double ry0 = node1->x[1] - node0->x[1];
        double rz0 = node1->x[2] - node0->x[2];
        double rx1 = node2->x[0] - node1->x[0];
        double ry1 = node2->x[1] - node1->x[1];
        double rz1 = node2->x[2] - node1->x[2];
        double nx = ry0 * rz1 - rz0 * ry1;
        double ny = rz0 * rx1 - rx0 * rz1;
        double nz = rx0 * ry1 - ry0 * rx1;
        double r = sqrt(nx * nx + ny * ny + nz * nz);

        element->surface_nx[iSurface][0] = nx / r;
        element->surface_nx[iSurface][1] = ny / r;
        element->surface_nx[iSurface][2] = nz / r;
      }
    }
  }
}
