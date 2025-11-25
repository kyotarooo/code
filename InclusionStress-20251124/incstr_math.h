#ifndef __INCSTR_MATH_H_
#define __INCSTR_MATH_H_

extern void INCSTR_NormalizeVector(double *a);

extern void INCSTR_NormalizeTensor(double **a);

extern double INCSTR_VectorLength(int n, double a[3]);

extern void INCSTR_VectorProduct(double a[3], double b[3], double c[3]);

extern void INCSTR_RotateCoordinateVector(double rotation_tensor[3][3],
                                          double coord[3],
                                          double rotated_coord[3]);

extern void INCSTR_RotateStressTensor(double rotation_tensor[3][3],
                                      double stress[3][3],
                                      double rotated_stress[3][3]);

extern void INCSTR_AddTensor(double a[3][3], double b[3][3]);

#endif