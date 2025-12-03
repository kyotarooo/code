#include <stdio.h>

#include "dd_element_handle.h"
#include "dd_node_handle.h"
#include "dd_struct.h"

static void RecordElementAnnihilation(int elementID0, int elementID1,
                                      DD_t* dd) {
  FILE* fp = dd->file.direct_interaction;
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

  fprintf(fp, "Time: %e\n", dd->step.t);
  fprintf(fp, "InteractionType: ElementAnnihilation\n");
  fprintf(fp, "ElementID0: %d NodeID0: %d type: %d NodeID1: %d type: %d\n",
          elementID0, nodeID00, node00->type, nodeID01, node01->type);
  fprintf(fp, "ElementID1: %d NodeID0: %d type: %d NodeID1: %d type: %d\n",
          elementID1, nodeID10, node10->type, nodeID11, node11->type);
  fflush(fp);
}

void DD_ElementAnnihilation(int elementID0, int elementID1, DD_t* dd) {
  ELEMENT_t* element0 = &dd->element[elementID0];
  ELEMENT_t* element1 = &dd->element[elementID1];
  int nodeID01 = element0->nodeID[1];
  int nodeID11 = element1->nodeID[1];

  RecordElementAnnihilation(elementID0, elementID1, dd);

  element0->nodeID[1] = nodeID11;
  element1->nodeID[1] = nodeID01;

  DD_TangentVector(elementID0, dd);
  DD_TangentVector(elementID1, dd);
}

static void RecordNodeAnnihilation(int nodeID0, int nodeID1, DD_t* dd) {
  FILE* fp = dd->file.direct_interaction;
  NODE_t* node0 = &dd->node[nodeID0];
  NODE_t* node1 = &dd->node[nodeID1];

  fprintf(fp, "Time: %e\n", dd->step.t);
  fprintf(fp, "InteractionType: NodeAnnihilation\n");
  fprintf(fp, "NodeID0: %d type: %d\n", nodeID0, node0->type);
  fprintf(fp, "NodeID1: %d type: %d\n", nodeID1, node1->type);
  fflush(fp);
}

void DD_NodeAnnihilation(int nodeID0, int nodeID1, DD_t* dd) {
  RecordNodeAnnihilation(nodeID0, nodeID1, dd);
  DD_ConnectElementPairs(nodeID0, nodeID1, dd);
  DD_ReplaceNode(nodeID1, nodeID0, dd);
  DD_DeleteNode(nodeID1, dd);

  DD_TangentVectors(dd);
}