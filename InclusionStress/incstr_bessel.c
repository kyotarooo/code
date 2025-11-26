#include <math.h>
#include <stdio.h>

double INCSTR_BesselFunction(int index1, int index2, int index3, double param[7]) {
  int index = 100 * (index1 + 1) + 10 * (index2 + 1) + 1 * (index3 + 1);
  double r = param[0];
  double sqr = param[1];
  double x = param[2];
  double m = param[3];
  double e = param[4];
  double f = param[5];
  double p = param[6];

  if (index == 220) {
    if (r < 1.0) {
      return (x / (2.0 * M_PI * sqr)) *
                 (2.0 / m * e - m / r * (1.0 + r * r + 0.5 * x * x) * f +
                  m / (2.0 * r) * (1.0 - r) * (1.0 - r) * p) +
             0.5 * r;
    } else {
      return (x / (2.0 * M_PI * sqr)) *
                 (2.0 / m * e - m / r * (1.0 + r * r + 0.5 * x * x) * f +
                  m / (2.0 * r) * (1.0 - r) * (1.0 - r) * p) +
             1.0 / (2.0 * r);
    }
  } else if (index == 210) {
    if (r < 1.0) {
      return (1.0 / (2.0 * M_PI * sqr)) *
                 (4.0 * r / m * e + m * (1.0 - r * r) * f +
                  m * x * x * (1.0 - r) / (1.0 + r) * p) -
             x;
    } else {
      return (1.0 / (2.0 * M_PI * sqr)) *
             (4.0 * r / m * e + m * (1.0 - r * r) * f +
              m * x * x * (1.0 - r) / (1.0 + r) * p);
    }
  } else if (index == 211) {
    if (r < 1.0) {
      return (-1.0 * m * x / (2.0 * M_PI * sqr)) *
                 (f + (1.0 - r) / (1.0 + r) * p) +
             1.0;
    } else {
      return (-1.0 * m * x / (2.0 * M_PI * sqr)) *
             (f + (1.0 - r) / (1.0 + r) * p);
    }
  } else if (index == 221) {
    return (1.0 / (M_PI * m * sqr)) * ((2.0 - m * m) * f - 2.0 * e);
  } else if (index == 231) {
    if (r < 1.0) {
      return (x / (M_PI * sqr * r)) *
             (2.0 / m * e - m / r * (1.0 + 0.5 * r * r + 0.5 * x * x) * f +
              m / (2.0 * r) * (1.0 - r) / (1.0 + r) * p);
    } else {
      return (x / (M_PI * sqr * r)) *
                 (2.0 / m * e - m / r * (1.0 + 0.5 * r * r + 0.5 * x * x) * f +
                  m / (2.0 * r) * (1.0 - r) / (1.0 + r) * p) +
             1.0 / (r * r);
    }
  }
  return 0.0;
}