#ifndef __DD_NODE_HANDLE_H_
#define __DD_NODE_HANDLE_H_

#include "dd_struct.h"

extern int DD_AddNode(int type, double x[3], double d[3], DD_t *dd);
extern void DD_DeleteNode(int nodeID, DD_t *dd);
extern void DD_ReplaceNode(int target_nodeID, int replace_nodeID, DD_t *dd);
extern void DD_FindElementsAroundNode(DD_t *dd);
extern double DD_ElementCosineAngle(int nodeID, DD_t *dd);
extern void DD_BackupPosition(DD_t *dd);

#endif
