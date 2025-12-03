#ifndef __DD_WRITE_H_
#define __DD_WRITE_H_

#include "dd_struct.h"

extern void DD_WriteLog(int stepID, DD_t *dd);
extern void DD_WriteTime(DD_t *dd, char *directory);
extern void DD_WriteDislocations(DD_t *dd, char *directory);
extern void DD_WriteMechanicalBehavior(DD_t *dd, char *directory);
extern void DD_WriteSimulationVolume(DD_t *dd, char *directory);

#endif
