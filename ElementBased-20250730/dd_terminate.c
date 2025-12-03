#include "dd_struct.h"
#include "dd_terminate.h"
#include <math.h>

int  DD_TerminateSimulation
(int step, DD_t *dd)
{
  double  goal_line = dd->bc.size[0]; // terminateの判断距離

  /***********************************/
  /** 全てのnodeが越えたらterminate **/
  /***********************************/

  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    double  x = dd->node[iNode].x[0];

    if (x < goal_line) return  0;
  }

  return  1;
}
