#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_defs.h"
#include "dd_node_handle.h"
#include "dd_periodic_bc.h"
#include "dd_struct.h"

static void SubstituteVector(double* a, double b[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] = a[iDim];
  }
}

static double ScalarProduct(double a[3], double b[3]) {
  double c = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    c += a[iDim] * b[iDim];
  }
  return c;
}

static void UnitVector(double a[3]) {
  double b = sqrt(ScalarProduct(a, a));

  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] /= b;
  }
}

void DD_ElementCenter(int elementID, DD_t* dd, double x[3]) {
  ELEMENT_t* element = &dd->element[elementID];
  int nodeID0 = element->nodeID[0];
  int nodeID1 = element->nodeID[1];
  NODE_t* node0 = &dd->node[nodeID0];
  NODE_t* node1 = &dd->node[nodeID1];
  double x0[3], x1[3];

  SubstituteVector(node0->x, x0);
  SubstituteVector(node1->x, x1);
  DD_PeriodicReturn(x0, x1, &dd->bc);

  for (int iDim = 0; iDim < 3; iDim++) {
    x[iDim] = (x0[iDim] + x1[iDim]) / 2.0;
  }
}

void DD_DeleteElement(int elementID, DD_t* dd) {
  int end_elementID = dd->nElements - 1;
  ELEMENT_t* element0 = &dd->element[elementID];
  ELEMENT_t* element1 = &dd->element[end_elementID];

  for (int iNode = 0; iNode < 2; iNode++) {
    element0->nodeID[iNode] = element1->nodeID[iNode];
  }
  for (int iDim = 0; iDim < 3; iDim++) {
    element0->burgers[iDim] = element1->burgers[iDim];
  }
  for (int iDim = 0; iDim < 3; iDim++) {
    element0->slip[iDim] = element1->slip[iDim];
  }
  element0->length = element1->length;
  for (int iDim = 0; iDim < 3; iDim++) {
    element0->tangent[iDim] = element1->tangent[iDim];
  }
  element0->coreID = element1->coreID;
  for (int iIntegral = 0; iIntegral < N_DD_INTEGRALS; iIntegral++) {
    INTEGRAL_t* integral0 = &element0->integral[iIntegral];
    INTEGRAL_t* integral1 = &element1->integral[iIntegral];

    for (int i = 0; i < dd->core.nDislocations; i++) {
      integral0->core_x[i] = integral1->core_x[i];
    }
  }

  dd->nElements -= 1;
}

int DD_TangentVector(int elementID, DD_t* dd) {
  ELEMENT_t* element = &dd->element[elementID];
  int nodeID0 = element->nodeID[0];
  int nodeID1 = element->nodeID[1];
  NODE_t* node0 = &dd->node[nodeID0];
  NODE_t* node1 = &dd->node[nodeID1];
  double x0[3], x1[3];
  double t[3];
  double length = 0.0;

  // Delete zero length element
  if (nodeID0 == nodeID1) {
    DD_DeleteElement(elementID, dd);
    if (nodeID0 != NO_NODE) {
      DD_DeleteNode(nodeID0, dd);
    }
    return 1;
  }

  for (int iDim = 0; iDim < 3; iDim++) {
    x0[iDim] = node0->x[iDim];
    x1[iDim] = node1->x[iDim];
  }
  DD_PeriodicReturn(x0, x1, &dd->bc);

  for (int iDim = 0; iDim < 3; iDim++) {
    t[iDim] = x1[iDim] - x0[iDim];

    length += t[iDim] * t[iDim];
  }
  element->length = length = sqrt(length);
  for (int iDim = 0; iDim < 3; iDim++) {
    element->tangent[iDim] = t[iDim] / length;
  }

  return 0;
}

void DD_TangentVectors(DD_t* dd) {
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    if (DD_TangentVector(iElement, dd)) {
      iElement -= 1;
    }
  }
}

int DD_AddElement(int nodeID[2], double burgers[3], double slip[3], int coreID,
                  DD_t* dd) {
  int elementID = dd->nElements;
  double r = 0.0;

  if (elementID == NMAX_ELEMENTS) {
    fprintf(stderr, "Number of elements exceeds over the maximum\n");
    exit(1);
  }

  for (int iNode = 0; iNode < 2; iNode++) {
    dd->element[elementID].nodeID[iNode] = nodeID[iNode];
  }
  for (int iDim = 0; iDim < 3; iDim++) {
    dd->element[elementID].burgers[iDim] = burgers[iDim];
  }
  r = sqrt(ScalarProduct(slip, slip));
  for (int iDim = 0; iDim < 3; iDim++) {
    dd->element[elementID].slip[iDim] = slip[iDim] / r;
  }
  dd->element[elementID].coreID = coreID;
  for (int iIntegral = 0; iIntegral < N_DD_INTEGRALS; iIntegral++) {
    INTEGRAL_t* integral = &dd->element[elementID].integral[iIntegral];

    for (int i = 0; i < dd->core.nDislocations; i++) {
      integral->core_x[i] = dd->core.property[coreID].x[i];
    }
  }
  dd->nElements += 1;
  DD_TangentVector(elementID, dd);

  return elementID;
}

void DD_DivideElement(int elementID, DD_t* dd) {
  ELEMENT_t* element = &dd->element[elementID];
  NODE_t* node0 = &dd->node[element->nodeID[0]];
  NODE_t* node1 = &dd->node[element->nodeID[1]];
  double x[3];
  double d[3];
  int new_nodeID;
  int nodeID[2];
  double burgers[3], slip[3];
  int coreID;
  int new_elementID;

  DD_ElementCenter(elementID, dd, x);
  SubstituteVector(node0->direction, d);

  // Add node at the middle of the element
  if (node0->type == node1->type) {
    new_nodeID = DD_AddNode(node0->type, x, d, dd);
  } else {
    new_nodeID = DD_AddNode(MOBILE_NODE, x, d, dd);
  }

  // Rearrange element
  nodeID[0] = new_nodeID;
  nodeID[1] = element->nodeID[1];
  element->nodeID[1] = new_nodeID;

  // Copy element data
  for (int iDim = 0; iDim < 3; iDim++) {
    burgers[iDim] = element->burgers[iDim];
    slip[iDim] = element->slip[iDim];
  }
  coreID = element->coreID;

  // Add element
  new_elementID = DD_AddElement(nodeID, burgers, slip, coreID, dd);

  // Tangent vectors
  DD_TangentVector(elementID, dd);
  DD_TangentVector(new_elementID, dd);
}

static int SameSlipSystem(ELEMENT_t* element0, ELEMENT_t* element1) {
  double b0[3], b1[3];
  double s0[3], s1[3];
  double burgers_cosine_angle;
  double slip_cosine_angle;
  double eps = 1.0e-06;

  SubstituteVector(element0->burgers, b0);
  SubstituteVector(element1->burgers, b1);
  UnitVector(b0);
  UnitVector(b1);
  SubstituteVector(element0->slip, s0);
  SubstituteVector(element1->slip, s1);

  burgers_cosine_angle = ScalarProduct(b0, b1);
  slip_cosine_angle = ScalarProduct(s0, s1);

  if (burgers_cosine_angle > 1.0 - eps && slip_cosine_angle > 1.0 - eps) {
    return 1;
  }
  return 0;
}

int DD_MobileElement(int elementID, DD_t* dd) {
  ELEMENT_t* element = &dd->element[elementID];
  int nMobile_nodes = 0;

  for (int iNode = 0; iNode < 2; iNode++) {
    int nodeID = element->nodeID[iNode];

    if (dd->node[nodeID].type == MOBILE_NODE) {
      nMobile_nodes += 1;
    }
  }

  if (nMobile_nodes == 0) {
    return 0;
  }
  return 1;
}

int ElementPairs(int nodeID0, int nodeID1, DD_t* dd,
                 int list[NMAX_ELEMENT_PAIRS][2]) {
  NODE_t* node0 = &dd->node[nodeID0];
  NODE_t* node1 = &dd->node[nodeID1];
  int nElements = node0->nElements;
  int nElement_pairs = 0;

  for (int iElement = 0; iElement < nElements; iElement++) {
    int elementID0 = node0->elementID[iElement];
    ELEMENT_t* element0 = &dd->element[elementID0];

    if (!DD_MobileElement(elementID0, dd)) continue;

    list[nElement_pairs][0] = elementID0;
    for (int jElement = 0; jElement < nElements; jElement++) {
      int elementID1 = node1->elementID[jElement];
      ELEMENT_t* element1 = &dd->element[elementID1];

      if (SameSlipSystem(element0, element1)) {
        list[nElement_pairs][1] = elementID1;
        break;
      }
    }
    nElement_pairs += 1;
  }
  return nElement_pairs;
}

void DD_ConnectElementPairs(int nodeID0, int nodeID1, DD_t* dd) {
  int element_pair[NMAX_ELEMENT_PAIRS][2];
  int nElement_pairs = ElementPairs(nodeID0, nodeID1, dd, element_pair);

  for (int iElement = 0; iElement < nElement_pairs; iElement++) {
    int elementID0 = element_pair[iElement][0];
    int elementID1 = element_pair[iElement][1];
    ELEMENT_t* element0 = &dd->element[elementID0];
    ELEMENT_t* element1 = &dd->element[elementID1];
    int nodeID00 = element0->nodeID[0];
    int nodeID01 = element0->nodeID[1];
    int nodeID10 = element1->nodeID[0];
    int nodeID11 = element1->nodeID[1];
    NODE_t* node00 = &dd->node[nodeID00];
    NODE_t* node01 = &dd->node[nodeID01];
    NODE_t* node10 = &dd->node[nodeID10];
    NODE_t* node11 = &dd->node[nodeID11];
    double x00[3], x01[3], x10[3], x11[3];
    double d[3] = {0.0, 0.0, 0.0};
    double x[3] = {0.0, 0.0, 0.0};
    int new_nodeID;

    SubstituteVector(node00->x, x00);
    SubstituteVector(node01->x, x01);
    SubstituteVector(node10->x, x10);
    SubstituteVector(node11->x, x11);

    DD_PeriodicReturn(x00, x01, &dd->bc);
    DD_PeriodicReturn(x00, x10, &dd->bc);
    DD_PeriodicReturn(x00, x11, &dd->bc);

    for (int iDim = 0; iDim < 3; iDim++) {
      x[iDim] = 0.25 * (x00[iDim] + x01[iDim] + x10[iDim] + x11[iDim]);
    }

    new_nodeID = DD_AddNode(MOBILE_NODE, x, d, dd);
    if (nodeID00 == nodeID0) {
      element0->nodeID[0] = new_nodeID;
    }
    if (nodeID01 == nodeID0) {
      element0->nodeID[1] = new_nodeID;
    }
    if (nodeID10 == nodeID1) {
      element1->nodeID[0] = new_nodeID;
    }
    if (nodeID11 == nodeID1) {
      element1->nodeID[1] = new_nodeID;
    }
  }
}

void DD_RemoveElement(int elementID, DD_t* dd) {
  ELEMENT_t* element = &dd->element[elementID];
  int target_nodeID = -1;
  int replace_nodeID = -1;
  int nodeID0 = element->nodeID[0];
  int nodeID1 = element->nodeID[1];
  NODE_t* node0 = &dd->node[nodeID0];
  NODE_t* node1 = &dd->node[nodeID1];

  // if (node0->type == MOBILENODE) {
  if (node0->nElements == 2) {
    target_nodeID = nodeID0;
    replace_nodeID = nodeID1;

    // if (node1->type == MOBILENODE) {
    if (node1->nElements == 2) {
      double x[3];

      DD_ElementCenter(elementID, dd, x);
      for (int iDim = 0; iDim < 3; iDim++) {
        node1->x[iDim] = x[iDim];
      }
    }
    DD_ReplaceNode(target_nodeID, replace_nodeID, dd);
    DD_DeleteNode(target_nodeID, dd);
    //} else if (node1->type == MOBILENODE) {
  } else if (node1->nElements == 2) {
    target_nodeID = nodeID1;
    replace_nodeID = nodeID0;

    // if (node0->type == MOBILENODE) {
    if (node0->nElements == 2) {
      double x[3];

      DD_ElementCenter(elementID, dd, x);
      for (int iDim = 0; iDim < 3; iDim++) {
        node0->x[iDim] = x[iDim];
      }
    }
    DD_ReplaceNode(target_nodeID, replace_nodeID, dd);
    DD_DeleteNode(target_nodeID, dd);
  } else if (node0->nElements == 1 && node1->nElements == 1) {
    // DD_DeleteElement(elementID, dd);
    DD_DeleteNode(nodeID0, dd);
    if (nodeID1 == dd->nNodes) nodeID1 = nodeID0;
    DD_DeleteNode(nodeID1, dd);
  } else {
    DD_ConnectElementPairs(nodeID0, nodeID1, dd);
    DD_DeleteNode(nodeID0, dd);
    if (nodeID1 == dd->nNodes) nodeID1 = nodeID0;
    DD_DeleteNode(nodeID1, dd);
  }

  // Replace the element
  DD_DeleteElement(elementID, dd);

  // Update tangent vector
  DD_TangentVectors(dd);
}

void DD_MergeElements(int elementID0, int elementID1, DD_t* dd) {
  ELEMENT_t* element0 = &dd->element[elementID0];
  ELEMENT_t* element1 = &dd->element[elementID1];
  int nodeID00 = element0->nodeID[0];
  int nodeID01 = element0->nodeID[1];
  int nodeID10 = element1->nodeID[0];
  int nodeID11 = element1->nodeID[1];
  double t0[3], t1[3];

  SubstituteVector(element0->tangent, t0);
  SubstituteVector(element1->tangent, t1);

  // Depends on the angle between the elements
  // Make new connections by replacing nodes
  if (ScalarProduct(t0, t1) < 0.0) {
    DD_ReplaceNode(nodeID11, nodeID00, dd);
    DD_ReplaceNode(nodeID10, nodeID01, dd);
  } else {
    DD_ReplaceNode(nodeID10, nodeID00, dd);
    DD_ReplaceNode(nodeID11, nodeID01, dd);
  }

  // Delete the element and the nodes
  DD_DeleteElement(elementID1, dd);
  // if (dd->node[nodeID10].type == MOBILE_NODE) {
  if (dd->node[nodeID10].nElements == 2) {
    DD_DeleteNode(nodeID10, dd);
    if (nodeID11 == dd->nNodes) nodeID11 = nodeID10;
  }
  // if (dd->node[nodeID11].type == MOBILE_NODE) {
  if (dd->node[nodeID11].nElements == 2) {
    DD_DeleteNode(nodeID11, dd);
  }
}