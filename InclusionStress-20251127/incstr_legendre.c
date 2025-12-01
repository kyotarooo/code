#include "incstr_legendre.h"

#include <math.h>
#include <stdio.h>

void INCSTR_LegendrePolynomials(double x, double p[NMAX_LEGENDRE + 1]) {
  p[0] = 1.0;
  p[1] = x;
  for (int i = 1; i < NMAX_LEGENDRE; i++) {
    double l = (double)i;

    p[i + 1] = 1.0 / (l + 1.0) * ((2.0 * l + 1.0) * x * p[i] - l * p[i - 1]);
  }
}

void INCSTR_AssociatedLegendrePolynomials(double x,
                                          double p1[NMAX_LEGENDRE + 1]) {
  double r = sqrt(1.0 - x * x);
  
  p1[0] = 0.0;
  p1[1] = -1.0 * r;
  for (int i = 1; i < NMAX_LEGENDRE; i++) {
    double l = (double)i;

    p1[i + 1] = 1.0 / l * ((2.0 * l + 1.0) * x * p1[i] - (l + 1.0) * p1[i - 1]);
  }
}

static double Bn(int n, double p[NMAX_LEGENDRE + 1]) {
  double c1 = 2.0 * (double)n - 1.0;
  double c2 = 2.0 * (double)n + 1.0;
  double c3 = 2.0 * (double)n - 3.0;

  if (n > 2) {
    return 1.0 / c1 *
           (1.0 / c2 * p[n + 1] - 2.0 * c1 / (c2 * c3) * p[n - 1] +
            1.0 / c3 * p[n - 3]);
  }
  return 1.0 / c1 * (1.0 / c2 * p[n + 1] - 2.0 * c1 / (c2 * c3) * p[n - 1]);
}

double INCSTR_LureBn(int n, double p1[NMAX_LEGENDRE + 1],
                     double p2[NMAX_LEGENDRE + 1]) {
  return Bn(n, p2) - Bn(n, p1);
}

static double Dn(int n, double p[NMAX_LEGENDRE + 1]) {
  double c1 = 2.0 * (double)n + 3.0;
  double c2 = 2.0 * (double)n + 5.0;
  double c3 = 2.0 * (double)n + 1.0;

  if (n > 0) {
    return 1.0 / c1 *
           (1.0 / c2 * p[n + 3] - 2.0 * c1 / (c2 * c3) * p[n + 1] +
            1.0 / c3 * p[n - 1]);
  }
  return 1.0 / c1 * (1.0 / c2 * p[n + 3] - 2.0 * c1 / (c2 * c3) * p[n + 1]);
}

double INCSTR_LureDn(int n, double p1[NMAX_LEGENDRE + 1],
                     double p2[NMAX_LEGENDRE + 1]) {
  return Dn(n, p2) - Dn(n, p1);
}