#ifndef __INCSTR_LEGENDRE_H_
#define __INCSTR_LEGENDRE_H_

#define NMAX_LEGENDRE 30

extern void INCSTR_LegendrePolynomials(double x, double p[NMAX_LEGENDRE + 1]);

extern void INCSTR_AssociatedLegendrePolynomials(double x,
                                                 double p1[NMAX_LEGENDRE + 1]);

extern double INCSTR_LureDn(int n, double p1[NMAX_LEGENDRE + 1],
                            double p2[NMAX_LEGENDRE + 1]);

extern double INCSTR_LureBn(int n, double p1[NMAX_LEGENDRE + 1],
                            double p2[NMAX_LEGENDRE + 1]);

#endif