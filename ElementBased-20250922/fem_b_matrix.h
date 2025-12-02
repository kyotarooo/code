#ifndef __FEM_B_MATRIX_H_
#define __FEM_B_MATRIX_H_

extern void FEM_MakeBMatrix(double b_matrix[6][60], double *jacobian,
                            double u[3], double element_x[20][3]);

#endif
