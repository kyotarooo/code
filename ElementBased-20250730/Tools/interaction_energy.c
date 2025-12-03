#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct __node_t_ {
  int type;
  double *x;
} NODE_t;

typedef struct __element_t_ {
  int *nodeID;
  double *burgers;
  double *slip;
  double *tangent;
} ELEMENT_t;

typedef struct __material_t_ {
  double g, v;
  double a_core;
} MATERIAL_t;

typedef struct __dd_t_ {
  int nNodes;
  NODE_t *node;
  int nElements;
  ELEMENT_t *element;
  MATERIAL_t material;

  double interaction_energy;
} DD_t;

static void UnitVector(double a[3]) {
  double r = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    r += a[iDim] * a[iDim];
  }
  r = sqrt(r);
  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] /= r;
  }
}

static void SubstituteVector(double *a, double b[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    b[iDim] = a[iDim];
  }
}

static void TangentVector(int elementID, DD_t *dd) {
  ELEMENT_t *element = &dd->element[elementID];
  int nodeID0 = element->nodeID[0];
  int nodeID1 = element->nodeID[1];
  NODE_t *node0 = &dd->node[nodeID0];
  NODE_t *node1 = &dd->node[nodeID1];
  double x0[3], x1[3];
  double t[3];

  SubstituteVector(node0->x, x0);
  SubstituteVector(node1->x, x1);
  for (int iDim = 0; iDim < 3; iDim++) {
    t[iDim] = x1[iDim] - x0[iDim];
  }
  UnitVector(t);
  for (int iDim = 0; iDim < 3; iDim++) {
    element->tangent[iDim] = t[iDim];
  }
}

static void ReadDislocations(DD_t *dd, char *file_name) {
  FILE *fp = fopen(file_name, "r");

  fscanf(fp, "%d", &dd->nNodes);
  dd->node = (NODE_t *)malloc(dd->nNodes * sizeof(NODE_t));
  for (int iNode = 0; iNode < dd->nNodes; iNode++) {
    NODE_t *node = &dd->node[iNode];

    node->x = (double *)malloc(3 * sizeof(double));
    fscanf(fp, "%d", &node->type);
    for (int iDim = 0; iDim < 3; iDim++) {
      fscanf(fp, "%lf", node->x[iDim]);
    }
  }

  fscanf(fp, "%d", &dd->nElements);
  dd->element = (ELEMENT_t *)malloc(dd->nElements * sizeof(ELEMENT_t));
  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    ELEMENT_t *element = &dd->element[iElement];
    double s[3];

    element->nodeID = (int *)malloc(2 * sizeof(int));
    element->burgers = (double *)malloc(3 * sizeof(double));
    element->slip = (double *)malloc(3 * sizeof(double));
    element->tangent = (double *)malloc(3 * sizeof(double));
    for (int iNode = 0; iNode < 2; iNode++) {
      fscanf(fp, "%d", &element->nodeID[iNode]);
    }
    for (int iDim = 0; iDim < 3; iDim++) {
      fscanf(fp, "%lf", &element->burgers[iDim]);
    }
    for (int iDim = 0; iDim < 3; iDim++) {
      fscanf(fp, "%lf", &element->slip[iDim]);
    }
    SubstituteVector(element->slip, s);
    UnitVector(s);
    for (int iDim = 0; iDim < 3; iDim++) {
      element->slip[iDim] = s[iDim];
    }
    TangentVector(iElement, dd);
  }
  fclose(fp);
}

static void VectorProduct(double a[3], double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static double ScalarProduct(double a[3], double b[3]) {
  double c = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    c += a[iDim] * b[iDim];
  }
  return c;
}

static double InteractionEnergySegments(int elementID0, int elementID1,
                                        DD_t *dd) {
  double e = 0.0;
  double g = dd->material.g;
  double p = dd->material.v;
  double w0 = g / (4.0 * M_PI * (1.0 - p));
  ELEMENT_t *element0 = &dd->element[elementID0];
  ELEMENT_t *element1 = &dd->element[elementID1];
  int nodeID00 = element0->nodeID[0];
  int nodeID01 = element0->nodeID[1];
  int nodeID10 = element1->nodeID[0];
  int nodeID11 = element1->nodeID[1];
  NODE_t *node00 = &dd->node[nodeID00];
  NODE_t *node01 = &dd->node[nodeID01];
  NODE_t *node10 = &dd->node[nodeID10];
  NODE_t *node11 = &dd->node[nodeID11];
  double x00[3], x01[3], x10[3], x11[3];
  double t0[3], t1[3];
  double b0[3], b1[3];
  double u[3], v[3], vp[3];
  double uu;
  double a[5], ap[5];
  double bt, bpt, btp, bptp, bbp, ttp, bu, bpu, bv, bpv, bvp, bpvp;

  SubstituteVector(node00->x, x00);
  SubstituteVector(node01->x, x01);
  SubstituteVector(node10->x, x10);
  SubstituteVector(node11->x, x11);
  SubstituteVector(element0->tangent, t0);
  SubstituteVector(element1->tangent, t1);
  SubstituteVector(element0->burgers, b0);
  SubstituteVector(element1->burgers, b1);

  VectorProduct(t0, t1, u);
  VectorProduct(u, t0, v);
  VectorProduct(t1, u, vp);

  uu = ScalarProduct(u, u);
  bt = ScalarProduct(b0, t0);
  bpt = ScalarProduct(b1, t0);
  btp = ScalarProduct(b0, t1);
  bptp = ScalarProduct(b1, t1);
  bbp = ScalarProduct(b0, b1);
  ttp = ScalarProdcut(t0, t1);
  bu = ScalarProduct(b0, u);
  bpu = ScalarProduct(b1, u);
  bv = ScalarProduct(b0, v);
  bpv = ScalarProduct(b1, v);
  bvp = ScalarProduct(b0, vp);
  bpvp = ScalarProduct(b1, vp);

  w0 /= uu;
  a[0] = (1.0 - p) * bt * bptp + 2.0 * p * bpt * btp;
  a[1] = (bbp + bt * bpt) * ttp;
  ap[1] = (bbp + btp * bptp) * ttp;
  ap[2] = (bu * bpvp + bpu * bvp) * ttp / uu;
  a[2] = (bu * bpv + bpu * bv) * ttp / uu;
  a[3] = (bt * bpv + btp * bpvp) * ttp;
  a[4] = 2.0 * (bu * bpu) * ttp / uu;
  
  return e;
}

static void InteractionEnergy(DD_t *dd) {
  dd->interaction_energy = 0.0;

  for (int iElement = 0; iElement < dd->nElements; iElement++) {
    for (int jElement = iElement + 1; jElement < dd->nElements; jElement++) {
      dd->interaction_energy +=
          InteractionEnergyElements(iElement, jElement, dd);
    }
  }
}
int main(int argc, char **argv) {
  DD_t dd;

  ReadDislocations(&dd, argv[1]);
  InteractionEnergy(&dd);

  return 0;
}