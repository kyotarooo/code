#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dd_defs.h"
#include "dd_element_handle.h"
#include "dd_node_handle.h"
#include "dd_periodic_bc.h"
#include "dd_struct.h"

static void SubstituteVector(double *a, double b[3]) {
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

static void OuterProduct(double a[3], double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static void UnitVector(double a[3]) {
  double r = sqrt(ScalarProduct(a, a));

  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] /= r;
  }
}

static void RelativePositionVector(double dx[3], double *x0, double *x1) {
  for (int iDim = 0; iDim < 3; iDim++) {
    dx[iDim] = x1[iDim] - x0[iDim];
  }
}

static void ControlDislocationElementLength(DD_t *dd) {
  double standard_element_length = dd->bc.element_length;

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    double length = dd->element[iElement].length;

    if (length > 2.0 * standard_element_length) {
      // If the element length is long, the element is divided into two
      DD_DivideElement(iElement, dd);
      DD_FindElementsAroundNode(dd);
      DD_TangentVectors(dd);
      iElement = -1;
      continue;
    } else if (length < 0.5 * standard_element_length) {
      // If the element length is short, the element is removed
      DD_RemoveElement(iElement, dd);
      DD_FindElementsAroundNode(dd);
      DD_TangentVectors(dd);
      iElement = -1;
      continue;
    }
  }
}

static void ControlDislocationElementAngle(DD_t *dd) {
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];
    double cosine_angle = DD_ElementCosineAngle(iNode, dd);

    if (cosine_angle >= 1.5) continue;

    // If the angle between elements around a node is small,
    // the node is removed
    if (cosine_angle > SMOOTH_COSINE_ANGLE) {
      int elementID1 = node->elementID[1];
      ELEMENT_t *element1 = &dd->element[elementID1];
      int target_nodeID = iNode;
      int replace_nodeID = element1->nodeID[0];

      if (replace_nodeID == iNode) {
        replace_nodeID = element1->nodeID[1];
      }

      DD_ReplaceNode(target_nodeID, replace_nodeID, dd);

      DD_DeleteNode(iNode, dd);
      DD_DeleteElement(elementID1, dd);

      DD_FindElementsAroundNode(dd);
      iNode = -1;
      continue;
    }
  }
}

#ifdef _FREE_SURFACE_
static int InsideVolume(int nodeID, DD_t *dd) {
  NODE_t *node = &dd->node[nodeID];
  double *size = dd->bc.size;
  double r = 0.01;
  double margin[3];

  for (int iDim = 0; iDim < 3; iDim++) {
    margin[iDim] = r * size[iDim];
  }

  for (int iDim = 0; iDim < 3; iDim++) {
    if (node->x[iDim] < -1.0 * margin[iDim] ||
        node->x[iDim] > size[iDim] + margin[iDim]) {
      return 0;
    }
  }
  return 1;
}

static void SetUpConstraintNodeOnFreeSurface(int nodeID, ELEMENT_t *element,
                                             DD_t *dd) {
  NODE_t *node = &dd->node[nodeID];
  double x[3];
  double s[3], n[3] = {0.0, 0.0, 0.0};
  double d[3];
  double size[3];
  double dmax = 0.0;
  int direction;

  SubstituteVector(node->x, x);
  SubstituteVector(element->slip, s);
  SubstituteVector(dd->bc.size, size);

  for (int iDim = 0; iDim < 3; iDim++) {
    double d = 0.0;

    if (x[iDim] < 0.0) {
      d = fabs(x[iDim]);
    } else if (x[iDim] > size[iDim]) {
      d = fabs(x[iDim] - size[iDim]);
    }
    if (dmax < d) {
      dmax = d;
      direction = iDim;
    }
  }

  n[direction] = 1.0;
  OuterProduct(s, n, d);
  UnitVector(d);
  for (int iDim = 0; iDim < 3; iDim++) {
    node->direction[iDim] = d[iDim];
  }
  // node->type = CONSTRAINT_NODE;
  node->type = FREE_SURFACE_NODE;
}

static void ArrangeNodesOnFreeSurface(DD_t *dd) {
  DD_FindElementsAroundNode(dd);

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];
    int nodeID0 = element->nodeID[0];
    int nodeID1 = element->nodeID[1];
    int outside[2] = {0, 0};
    for (int iNode = 0; iNode < 2; iNode++) {
      int nodeID = element->nodeID[iNode];

      if (!InsideVolume(nodeID, dd)) {
        SetUpConstraintNodeOnFreeSurface(nodeID, element, dd);
        outside[iNode] = 1;
      }
    }

    // If the element is outside, the element is deleted.
    if (outside[0] + outside[1] == 2) {
      DD_DeleteElement(iElement, dd);
      // It the node is edge node, the node must be deleted.
      if (dd->node[nodeID0].nElements == 1) {
        DD_DeleteNode(nodeID0, dd);
        if (nodeID1 > nodeID0) {
          nodeID1 = nodeID1 - 1;
        }
      }
      if (dd->node[nodeID1].nElements == 1) {
        DD_DeleteNode(nodeID1, dd);
      }
      DD_FindElementsAroundNode(dd);
      iElement = -1;
    } else if (outside[0]) {
      if (dd->node[nodeID0].nElements == 1 &&
          dd->node[nodeID1].type != MOBILE_NODE) {
        DD_DeleteElement(iElement, dd);
        DD_DeleteNode(nodeID0, dd);
        DD_FindElementsAroundNode(dd);
        iElement = -1;
      }
    } else if (outside[1]) {
      if (dd->node[nodeID1].nElements == 1 &&
          dd->node[nodeID0].type != MOBILE_NODE) {
        DD_DeleteElement(iElement, dd);
        DD_DeleteNode(nodeID1, dd);
        DD_FindElementsAroundNode(dd);
        iElement = -1;
      }
    }
  }
}

static double CosineAngleFreeSurface(double t[3], NODE_t *node) {
  double d[3];

  SubstituteVector(node->direction, d);
  return ScalarProduct(t, d);
}

static void ControlDislocationElementFreeSurfaceAngle(DD_t *dd) {
  double scale = 0.5;

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];
    double t[3];

    SubstituteVector(element->tangent, t);
    for (int iNode = 0; iNode < 2; iNode++) {
      int nodeID = element->nodeID[iNode];
      NODE_t *node = &dd->node[nodeID];
      int type = node->type;

      if (type == FREE_SURFACE_NODE) {
        double cos_angle = CosineAngleFreeSurface(t, node);

        // If the angle between the element and the free surface is small
        // enough, the element should be absorved by the free surface. To
        // simulate the behavior, the free surface node is slightly moved along
        // the free surface for reducing the angle.
        if (fabs(cos_angle) > FREE_SURFACE_COSINE_ANGLE) {
          int nodeID0 = element->nodeID[0];
          int nodeID1 = element->nodeID[1];
          NODE_t *node0 = &dd->node[nodeID0];
          NODE_t *node1 = &dd->node[nodeID1];

          if (node0->type == FREE_SURFACE_NODE &&
              node1->type == FREE_SURFACE_NODE) {
            DD_DeleteNode(nodeID0, dd);
            if (nodeID1 > nodeID0) nodeID1 -= 1;
            DD_DeleteNode(nodeID1, dd);
            DD_DeleteElement(iElement, dd);
          } else if (node0->type == FREE_SURFACE_NODE) {
            double dx[3], d[3], u;

            SubstituteVector(node0->direction, d);
            RelativePositionVector(dx, node0->x, node1->x);
            u = ScalarProduct(d, dx);
            for (int iDim = 0; iDim < 3; iDim++) {
              node0->x[iDim] += scale * u * d[iDim];
            }
          } else {
            double dx[3], d[3], u;

            SubstituteVector(node1->direction, d);
            RelativePositionVector(dx, node1->x, node0->x);
            u = ScalarProduct(d, dx);
            for (int iDim = 0; iDim < 3; iDim++) {
              node1->x[iDim] += scale * u * d[iDim];
            }
          }
        }
        break;
      }
    }
  }
}
#endif

void DD_RearrangeDislocations(DD_t *dd) {
  DD_TangentVectors(dd);
  DD_FindElementsAroundNode(dd);

  ControlDislocationElementAngle(dd);
  ControlDislocationElementLength(dd);

  DD_PeriodicBoundaryCondition(dd);

#ifdef _FREE_SURFACE_
  ArrangeNodesOnFreeSurface(dd);
  ControlDislocationElementFreeSurfaceAngle(dd);
#endif
}
