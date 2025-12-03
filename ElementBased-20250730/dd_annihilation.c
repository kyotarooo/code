#include <stdio.h>

#include "dd_element_handle.h"
#include "dd_node_handle.h"
#include "dd_struct.h"

void DD_ElementAnnihilation(int elementID0, int elementID1, DD_t *dd) {
  ELEMENT_t *element0 = &dd->element[elementID0];
  ELEMENT_t *element1 = &dd->element[elementID1];
  int nodeID01 = element0->nodeID[1];
  int nodeID11 = element1->nodeID[1];

  element0->nodeID[1] = nodeID11;
  element1->nodeID[1] = nodeID01;

  DD_TangentVector(elementID0, dd);
  DD_TangentVector(elementID1, dd);
}

void DD_NodeAnnihilation(int nodeID0, int nodeID1, DD_t *dd) {
  DD_ConnectElementPairs(nodeID0, nodeID1, dd);
  DD_ReplaceNode(nodeID1, nodeID0, dd);
  DD_DeleteNode(nodeID1, dd);

  DD_TangentVectors(dd);
}