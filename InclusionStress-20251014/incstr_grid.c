#include <stdio.h>

int INCSTR_GridID(int x[3], int *n) { return n[2] * (n[1] * x[0] + x[1]) + x[2]; }