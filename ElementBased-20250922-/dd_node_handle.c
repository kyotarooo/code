#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_defs.h"
#include "dd_periodic_bc.h"
#include "dd_struct.h"

static void SubstituteVector(double *a, double b[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] = a[iDim];
  }
}

static void OppositeVector(double a[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] *= -1.0;
  }
}

static double ScalarProduct(double a[3], double b[3]) {
  double c = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    c += a[iDim] * b[iDim];
  }
  return c;
}

int DD_AddNode(int type, double x[3], double d[3], DD_t *dd) {
  int nodeID = dd->nNodes;

  if (nodeID == NMAX_NODES) {
    fprintf(stderr, "Number of nodes exceeds over the maximum.\n");
    exit(1);
  }

  dd->node[nodeID].type = type;
  dd->node[nodeID].nElements = 0;
  for (int iDim = 0; iDim < 3; iDim++) {
    dd->node[nodeID].x[iDim] = x[iDim];
    dd->node[nodeID].direction[iDim] = d[iDim];
  }

  dd->nNodes += 1;
  return nodeID;
}

void DD_DeleteNode(int nodeID, DD_t *dd) {
  for (int iNode = nodeID; iNode < dd->nNodes - 1; iNode++) {
    NODE_t *node0 = &dd->node[iNode];
    NODE_t *node1 = &dd->node[iNode + 1];

    node0->type = node1->type;
    node0->nElements = node1->nElements;
    for (int iElement = 0; iElement < node0->nElements; iElement++) {
      node0->elementID[iElement] = node1->elementID[iElement];
    }
    for (int iDim = 0; iDim < 3; iDim++) {
      node0->x[iDim] = node1->x[iDim];
      node0->direction[iDim] = node1->direction[iDim];
    }
  }

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    for (int iNode = 0; iNode < 2; iNode++) {
      int element_nodeID = element->nodeID[iNode];

      if (element_nodeID > nodeID) {
        element->nodeID[iNode] -= 1;
      } else if (element_nodeID == nodeID) {
        element->nodeID[iNode] = NO_NODE;
      }
    }
  }

  dd->nNodes -= 1;
}

void DD_ReplaceNode(int target_nodeID, int replace_nodeID, DD_t *dd) {
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    for (int iNode = 0; iNode < 2; iNode++) {
      int nodeID = element->nodeID[iNode];

      if (nodeID == target_nodeID) {
        element->nodeID[iNode] = replace_nodeID;
      }
    }
  }
}

void DD_FindElementsAroundNode(DD_t *dd) {
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    dd->node[iNode].nElements = 0;
  }
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];

    for (int iNode = 0; iNode < 2; iNode++) {
      int nodeID = element->nodeID[iNode];
      NODE_t *node = &dd->node[nodeID];
      int elementID = node->nElements;

      if (elementID == NMAX_ELEMENTS_AROUND_NODE) {
        fprintf(stderr,
                "Number of elements around a node exceeds over the maximum.\n");
        exit(1);
      }
      node->elementID[elementID] = iElement;
      node->nElements += 1;
    }
  }

  // Delete isolated node
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];

    if (node->nElements == 0) {
      DD_DeleteNode(iNode, dd);
      iNode -= 1;
    }
  }
}

double DD_ElementCosineAngle(int nodeID, DD_t *dd) {
  NODE_t *node = &dd->node[nodeID];
  int elementID0, elementID1;
  ELEMENT_t *element0, *element1;
  double t0[3], t1[3];

  if (node->nElements != 2) return 2.0;

  elementID0 = node->elementID[0];
  elementID1 = node->elementID[1];
  element0 = &dd->element[elementID0];
  element1 = &dd->element[elementID1];

  // Tangent vector of elementID0
  SubstituteVector(element0->tangent, t0);
  if (element0->nodeID[1] == nodeID) {
    OppositeVector(t0);
  }

  // Tangent vector of elementID1
  SubstituteVector(element1->tangent, t1);
  if (element1->nodeID[1] == nodeID) {
    OppositeVector(t1);
  }

  // Calculate cosine angle between the elements
  return ScalarProduct(t0, t1);
}

void DD_BackupPosition(DD_t *dd) {
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];
    for (int iDim = 0; iDim < 3; iDim++) {
      node->x0[iDim] = node->x[iDim];
    }
  }
}