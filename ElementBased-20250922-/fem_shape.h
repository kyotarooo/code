#ifndef __FEM_SHAPE_H_
#define __FEM_SHAPE_H_

extern void FEM_ShapeFunction3D(double u[3], double shape[20]);
extern void FEM_DShapeFunction3D(double u[3], double dshape[20][3]);
extern void FEM_ShapeFunction2D(double u[2], double shape[8]);
extern void FEM_DShapeFunction2D(double u[2], double dshape[8][2]);
#endif
