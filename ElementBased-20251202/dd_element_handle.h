#ifndef __DD_ELEMENT_HANDLE_H_
#define __DD_ELEMENT_HANDLE_H_

#include "dd_defs.h"
#include "dd_struct.h"

extern void DD_ElementCenter(int elementID, DD_t *dd, double x[3]);
extern int DD_TangentVector(int elementID, DD_t *dd);
extern void DD_TangentVectors(DD_t *dd);
extern int DD_AddElement(int nodeID[2], double burgers[3], double slip[3],
                         int coreID, DD_t *dd);
extern void DD_DivideElement(int elementID, DD_t *dd);
extern void DD_DeleteElement(int elementID, DD_t *dd);
extern void DD_ConnectElementPairs(int nodeID0, int nodeID1, DD_t *dd);
extern int DD_MobileElement(int elementID, DD_t *dd);
extern void DD_RemoveElement(int elementID, DD_t *dd);
extern void DD_MergeElements(int elementID0, int elementID1, DD_t *dd);

#endif
