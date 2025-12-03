#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_annihilation.h"
#include "dd_defs.h"
#include "dd_element_handle.h"
#include "dd_junction.h"
#include "dd_node_handle.h"
#include "dd_periodic_bc.h"
#include "dd_struct.h"

static double ScalarProduct(double a[3], double b[3]) {
  double c = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    c += a[iDim] * b[iDim];
  }
  return c;
}

static void SubstituteVector(double *a, double b[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] = a[iDim];
  }
}

static void RelativeVector(double a[3], double b[3], double c[3]) {
  for(int iDim = 0; iDim < 3; iDim++) {
    c[iDim] = b[iDim] - a[iDim];
  }
}

static double SquaredDistanceBetweenElements(int elementID0, int elementID1,
                                             DD_t *dd) {
  ELEMENT_t *element0 = &dd->element[elementID0];
  double x0[3], x1[3];
  double t[3], R[3];
  double d[3];
  double Rdt;

  DD_ElementCenter(elementID0, dd, x0);
  DD_ElementCenter(elementID1, dd, x1);
  DD_PeriodicReturn(x0, x1, &dd->bc);

  RelativeVector(x0, x1, R);
  SubstituteVector(element0->tangent, t);

  Rdt = ScalarProduct(R, t);
  if (fabs(Rdt) > 0.5 * element0->length) return 1.0e+20;

  for (int iDim = 0; iDim < 3; iDim++) {
    d[iDim] = R[iDim] - Rdt * t[iDim];
  }
  return ScalarProduct(d, d);
}

static int NeighborElement(int elementID0, int elementID1, DD_t *dd) {
  ELEMENT_t *element0 = &dd->element[elementID0];
  ELEMENT_t *element1 = &dd->element[elementID1];

  for (int iNode = 0; iNode < 2; iNode++) {
    int nodeID0 = element0->nodeID[iNode];

    // if (dd->node[nodeID0].type != MOBILE_NODE) continue;
    if (dd->node[nodeID0].nElements != 2) continue;
    for (int jNode = 0; jNode < 2; jNode++) {
      int nodeID1 = element1->nodeID[jNode];

      // if (dd->node[nodeID1].type != MOBILE_NODE) continue;
      if (dd->node[nodeID1].nElements != 2) continue;
      if (nodeID0 == nodeID1) {
        return 1;
      }
    }
  }

  return 0;
}


static int FindNearestElement(int elementID, DD_t *dd, double *distance) {
  int nearest_elementID = NO_ELEMENT;
  double d2min = 1.0e+20;

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    double d2;

    if (iElement == elementID) continue;
    if (NeighborElement(elementID, iElement, dd)) continue;
    if (!DD_MobileElement(iElement, dd)) continue;

    d2 = SquaredDistanceBetweenElements(elementID, iElement, dd);
    if (d2min > d2) {
      d2min = d2;
      nearest_elementID = iElement;
    }
  }
  (*distance) = sqrt(d2min);
  return nearest_elementID;
}

static int SlipPlanePair(double *slip0, double *slip1) {
  double scalar = 0.0;
  double eps = 0.99;

  for (int iDim = 0; iDim < 3; iDim++) {
    scalar += slip0[iDim] * slip1[iDim];
  }

  if (fabs(scalar) > eps) return 1;

  return 0;
}

static int AttractiveInteraction(double b0[3], double t0[3], double b1[3],
                                 double t1[3]) {
  double b02 = ScalarProduct(b0, b0);
  double b12 = ScalarProduct(b1, b1);
  double b2[3];
  double b22;

  // Apply Frank rule to check if the interaction is attractive
  if (ScalarProduct(t0, t1) < 0.0) {
    for (int iDim = 0; iDim < 3; iDim++) {
      b2[iDim] = b0[iDim] - b1[iDim];
    }
  } else {
    for (int iDim = 0; iDim < 3; iDim++) {
      b2[iDim] = b0[iDim] + b1[iDim];
    }
  }
  b22 = ScalarProduct(b2, b2);

  if (b22 < b02 + b12) {
    return 1;
  }
  return 0;
}

static int InteractionType(int elementID0, int elementID1, DD_t *dd) {
  ELEMENT_t *element0 = &dd->element[elementID0];
  ELEMENT_t *element1 = &dd->element[elementID1];
  int same_slip_plane = SlipPlanePair(element0->slip, element1->slip);
  double b0[3], b1[3];
  double t0[3], t1[3];

  SubstituteVector(element0->burgers, b0);
  SubstituteVector(element1->burgers, b1);
  SubstituteVector(element0->tangent, t0);
  SubstituteVector(element1->tangent, t1);

  if (fabs(ScalarProduct(t0, t1)) < DIRECT_INTERACTION_COSINE_ANGLE) {
    return NO_INTERACTION;
  }
  if (!AttractiveInteraction(b0, t0, b1, t1)) {
    return NO_INTERACTION;
  }

  if (same_slip_plane) {
    return ANNIHILATION;
  } else {
    return JUNCTION;
  }

  return NO_INTERACTION;
}

static void ElementBased(DD_t *dd) {
  double rc = dd->bc.r_interaction;

  DD_TangentVectors(dd);
  DD_FindElementsAroundNode(dd);
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    double distance;
    int nearest_elementID;
    int type;

    if (!DD_MobileElement(iElement, dd)) continue;

    nearest_elementID = FindNearestElement(iElement, dd, &distance);
    if (nearest_elementID == NO_ELEMENT) continue;
    if (distance > rc) continue;

    type = InteractionType(iElement, nearest_elementID, dd);
    switch (type) {
      case ANNIHILATION:
        DD_ElementAnnihilation(iElement, nearest_elementID, dd);
        break;
      case JUNCTION:
        DD_Junction(iElement, nearest_elementID, dd);
        break;
    }
  }
}

static double SquaredDistanceBetweenNodes(int nodeID0, int nodeID1, DD_t *dd) {
  double x0[3], x1[3];
  double d2 = 0.0;

  SubstituteVector(dd->node[nodeID0].x, x0);
  SubstituteVector(dd->node[nodeID1].x, x1);
  DD_PeriodicReturn(x0, x1, &dd->bc);

  for (int iDim = 0; iDim < 3; iDim++) {
    double dx = x1[iDim] - x0[iDim];

    d2 += dx * dx;
  }
  return d2;
}

static int FindNearestNode(int nodeID, DD_t *dd, double *distance) {
  int nearest_nodeID = NO_NODE;
  double d2min = 1.0e+20;

  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];

    if (iNode == nodeID) continue;
    if (node->nElements != 2) {
      int mobile_element = 0;
      for (int iElement = 0; iElement < node->nElements; iElement++) {
        if (DD_MobileElement(node->elementID[iElement], dd)) {
          mobile_element = 1;
          break;
        }
      }
      if (mobile_element) {
        double d2 = SquaredDistanceBetweenNodes(nodeID, iNode, dd);

        if (d2min > d2) {
          d2min = d2;
          nearest_nodeID = iNode;
        }
      }
    }
  }
  (*distance) = sqrt(d2min);
  return nearest_nodeID;
}

static void NodeBased(DD_t *dd) {
  double rc = dd->bc.r_interaction;

  DD_FindElementsAroundNode(dd);
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];

    if (node->nElements != 2) {
      int mobile_element = 0;
      int nearest_nodeID;
      double distance;

      for (int iElement = 0; iElement < node->nElements; iElement++) {
        if (DD_MobileElement(node->elementID[iElement], dd)) {
          mobile_element = 1;
          break;
        }
      }

      if (!mobile_element) continue;

      nearest_nodeID = FindNearestNode(iNode, dd, &distance);
      if (nearest_nodeID == NO_NODE) continue;
      if (distance > rc) continue;

      DD_NodeAnnihilation(iNode, nearest_nodeID, dd);
      DD_FindElementsAroundNode(dd);
      iNode = -1;
    }
  }
}

void DD_DirectInteractions(DD_t *dd) {
  ElementBased(dd);
  NodeBased(dd);
}
