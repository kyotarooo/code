
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct _section_t_ {
  double *x;      // 断面上の一点の座標
  double *n;      // 断面の法線ベクトル
  double *s, *m;  // 応力成分（すべり系）
} SECTION_t;

typedef struct __grid_t_ {
  int *nx;           // 各方向のグリッド数
  double *size;      // 各方向のグリッドの大きさ
  double ***stress;  // グリッド点における応力
} GRID_t;

typedef struct __volume_t_ {
  SECTION_t section;
  GRID_t grid;
  double *size;  // 各方向のボリュームの大きさ
  int *nx;       // 各方向のグリッドセル数
} VOLUME_t;

static void UnitVector(double *a) {
  double r = 0.0;

  for (int iDim = 0; iDim < 3; iDim++) {
    r += a[iDim] * a[iDim];
  }
  r = sqrt(r);
  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] /= r;
  }
}

static double StressComponent(SECTION_t *section, double local_stress[9]) {
  int index[3][3] = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};
  double v = 0.0;
  double *s = section->s;
  double *m = section->m;

  // 入力された力の方向の面の方向から応力成分を計算する
  for (int iDim = 0; iDim < 3; iDim++) {
    for (int jDim = 0; jDim < 3; jDim++) {
      v += local_stress[index[iDim][jDim]] * s[iDim] * m[jDim];
    }
  }
  return v;
}

static void ReadVolume(VOLUME_t *volume, char *directory) {
  FILE *fp;
  char file_name[256];

  sprintf(file_name, "%s/section.inp", directory);
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot open %s file\n", file_name);
    exit(1);
  }

  // 領域のサイズ
  volume->size = (double *)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%lf", &volume->size[iDim]);
  }

  // 断面上の一点の座標
  volume->section.x = (double *)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%lf", &volume->section.x[iDim]);
  }

  // 断面の法線ベクトル
  volume->section.n = (double *)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%lf", &volume->section.n[iDim]);
  }
  UnitVector(volume->section.n);

  // すべり方向
  volume->section.s = (double *)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%lf", &volume->section.s[iDim]);
  }
  UnitVector(volume->section.s);

  // すべり面方向
  volume->section.m = (double *)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%lf", &volume->section.m[iDim]);
  }
  UnitVector(volume->section.m);

  // グリッドセルの数
  volume->nx = (int *)malloc(3 * sizeof(int));
  for (int iDim = 0; iDim < 3; iDim++) {
    fscanf(fp, "%d", &volume->nx[iDim]);
  }
  fclose(fp);
}

static void ReadSubgrids(int depth, int nx[100][3], FILE *fp) {
  for(int ix = 0; ix < nx[depth][0]; ix++) {
    for(int iy = 0; iy < nx[depth][1]; iy++) {
      for(int iz = 0; iz < nx[depth][2]; iz++) {
        double stress[9];
        int have_subgrid = 0;

        fread(stress, sizeof(double), 9, fp);
        fread(&have_subgrid, sizeof(int), 1, fp);
        if(have_subgrid) {
          ReadSubgrids(depth + 1, nx, fp);
        }
      }
    }
  }
}

static void ReadStress(VOLUME_t *volume, char *directory) {
  FILE *fp;
  char file_name[256];
  int *nx = volume->nx;
  double *size = volume->size;
  GRID_t *grid = &volume->grid;
  SECTION_t *section = &volume->section;
  int max_depth;
  int nx_depth[100][3];

  sprintf(file_name, "%s/grid_stress.inp", directory);
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot open %s file\n", file_name);
    exit(1);
  }

  // サブグリッドの深さ
  fread(&max_depth, sizeof(int), 1, fp);

  // 各方向のグリッド点の数
  for(int i = 0; i < max_depth; i++) {
    fread(nx_depth[i], sizeof(int), 3, fp);
  }
  grid->nx = (int *)malloc(3 * sizeof(int));
  for(int iDim = 0; iDim < 3; iDim++) {
    grid->nx[iDim] = nx_depth[0][iDim];
  }

  // 各方向のグリッドの大きさ
  grid->size = (double *)malloc(3 * sizeof(double));
  for (int iDim = 0; iDim < 3; iDim++) {
    grid->size[iDim] =
        size[iDim] / (double)nx[iDim] / ((double)grid->nx[iDim] - 1.0);
  }

  // 応力の読み込み
  grid->stress = (double ***)malloc(grid->nx[0] * sizeof(double **));
  for (int i = 0; i < grid->nx[0]; i++) {
    grid->stress[i] = (double **)malloc(grid->nx[1] * sizeof(double *));
    for (int j = 0; j < grid->nx[1]; j++) {
      grid->stress[i][j] = (double *)malloc(grid->nx[2] * sizeof(double));
      for (int k = 0; k < grid->nx[2]; k++) {
        double local_stress[9];
        int have_subgrid;

        fread(local_stress, sizeof(double), 9, fp);
        grid->stress[i][j][k] = StressComponent(section, local_stress);

        fread(&have_subgrid, sizeof(int), 1, fp);
        if(have_subgrid) {
          ReadSubgrids(1, nx_depth, fp);
        }
      }
    }
  }
  fclose(fp);
}

static void VTKHeader(FILE *fp) {
  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "Stress distribution on section\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
}

static double FindX(SECTION_t *s, double yi, double zi) {
  return (s->n[1] * (s->x[1] - yi) + s->n[2] * (s->x[2] - zi)) / s->n[0] +
         s->x[0];
}

static double FindY(SECTION_t *s, double xi, double zi) {
  return (s->n[0] * (s->x[0] - xi) + s->n[2] * (s->x[2] - zi)) / s->n[1] +
         s->x[1];
}

static double FindZ(SECTION_t *s, double xi, double yi) {
  return (s->n[0] * (s->x[0] - xi) + s->n[1] * (s->x[1] - yi)) / s->n[2] +
         s->x[2];
}

static void ShapeFunction(double u[3], double shape[8]) {
  double ux = u[0];
  double uy = u[1];
  double uz = u[2];

  shape[0] = 0.125 * (1.0 - ux) * (1.0 - uy) * (1.0 - uz);
  shape[1] = 0.125 * (1.0 + ux) * (1.0 - uy) * (1.0 - uz);
  shape[2] = 0.125 * (1.0 + ux) * (1.0 + uy) * (1.0 - uz);
  shape[3] = 0.125 * (1.0 - ux) * (1.0 + uy) * (1.0 - uz);
  shape[4] = 0.125 * (1.0 - ux) * (1.0 - uy) * (1.0 + uz);
  shape[5] = 0.125 * (1.0 + ux) * (1.0 - uy) * (1.0 + uz);
  shape[6] = 0.125 * (1.0 + ux) * (1.0 + uy) * (1.0 + uz);
  shape[7] = 0.125 * (1.0 - ux) * (1.0 + uy) * (1.0 + uz);
}

static void Offset(VOLUME_t *volume, double x[3]) {
  for (int iDim = 0; iDim < 3; iDim++) {
    double size = volume->size[iDim] / (double)volume->nx[iDim];
    int id = (int)floor(x[iDim] / size);

    if (id < 0)
      id = 0;
    else if (id >= volume->nx[iDim])
      id = volume->nx[iDim] - 1;

    x[iDim] -= (double)id * size;
  }
}

static double Stress(VOLUME_t *volume, double x[3]) {
  GRID_t *grid = &volume->grid;
  int gridID[3];
  double u[3];
  double shape[8];
  int grid_index[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                          {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  double local_stress = 0.0;

  // 対象となるグリッドセルの位置に合わせる
  Offset(volume, x);

  // 対象となるグリッドの検索とグリッド内座標の決定
  for (int iDim = 0; iDim < 3; iDim++) {
    gridID[iDim] = (int)floor(x[iDim] / grid->size[iDim]);
    if (gridID[iDim] < 0)
      gridID[iDim] = 0;
    else if (gridID[iDim] >= grid->nx[iDim] - 1) {
      gridID[iDim] = grid->nx[iDim] - 2;
    }

    u[iDim] = 2.0 * (x[iDim] / grid->size[iDim] - (double)gridID[iDim]) - 1.0;
    if (u[iDim] < -1.0)
      u[iDim] = -1.0;
    else if (u[iDim] > 1.0)
      u[iDim] = 1.0;
  }

  // 形状関数（六面体一次要素）を用いて応力を内挿する
  ShapeFunction(u, shape);
  for (int i = 0; i < 8; i++) {
    double s = shape[i];
    int target_gridID[3];

    for (int iDim = 0; iDim < 3; iDim++) {
      target_gridID[iDim] = gridID[iDim] + grid_index[i][iDim];
    }
    local_stress +=
        s * grid->stress[target_gridID[0]][target_gridID[1]][target_gridID[2]];
  }
  return local_stress;
}

static void XYSection(VOLUME_t *volume, FILE *fp) {
  SECTION_t *section = &volume->section;
  GRID_t *grid = &volume->grid;
  int nx[3] = {volume->nx[0] * (grid->nx[0] - 1) + 1,
               volume->nx[1] * (grid->nx[1] - 1) + 1,
               volume->nx[2] * (grid->nx[2] - 1) + 1};
  double *size = grid->size;
  int nPoints = nx[0] * nx[1];
  int nCells = (nx[0] - 1) * (nx[1] - 1);

  // 点の座標の書き出し
  fprintf(fp, "POINTS %d float\n", nPoints);
  for (int ix = 0; ix < nx[0]; ix++) {
    double xi = (double)ix * size[0];
    for (int iy = 0; iy < nx[1]; iy++) {
      double yi = (double)iy * size[1];
      double zi = FindZ(section, xi, yi);

      fprintf(fp, "%e %e %e\n", xi, yi, zi);
    }
  }

  // セルを構成する点番号の書き出し
  fprintf(fp, "CELLS %d %d\n", nCells, 5 * nCells);
  for (int ix = 0; ix < nx[0] - 1; ix++) {
    for (int iy = 0; iy < nx[1] - 1; iy++) {
      int nodeID0 = nx[1] * (ix + 0) + iy + 0;
      int nodeID1 = nx[1] * (ix + 1) + iy + 0;
      int nodeID2 = nx[1] * (ix + 1) + iy + 1;
      int nodeID3 = nx[1] * (ix + 0) + iy + 1;

      fprintf(fp, "4 %d %d %d %d\n", nodeID0, nodeID1, nodeID2, nodeID3);
    }
  }

  // セルタイプの書き出し（9: Rectangular）
  fprintf(fp, "CELL_TYPES %d\n", nCells);
  for (int i = 0; i < nCells; i++) {
    fprintf(fp, "9\n");
  }

  // 応力の書き出し
  fprintf(fp, "POINT_DATA %d\n", nPoints);
  fprintf(fp, "SCALARS Stress float\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (int ix = 0; ix < nx[0]; ix++) {
    double xi = (double)ix * size[0];
    for (int iy = 0; iy < nx[1]; iy++) {
      double yi = (double)iy * size[1];
      double zi = FindZ(section, xi, yi);
      double x[3] = {xi, yi, zi};

      fprintf(fp, "%e\n", Stress(volume, x));
    }
  }
}

static void YZSection(VOLUME_t *volume, FILE *fp) {
  SECTION_t *section = &volume->section;
  GRID_t *grid = &volume->grid;
  int nx[3] = {volume->nx[0] * (grid->nx[0] - 1) + 1,
               volume->nx[1] * (grid->nx[1] - 1) + 1,
               volume->nx[2] * (grid->nx[2] - 1) + 1};
  double *size = grid->size;
  int nPoints = nx[1] * nx[2];
  int nCells = (nx[1] - 1) * (nx[2] - 1);

  // 点の座標の書き出し
  fprintf(fp, "POINTS %d float\n", nPoints);
  for (int iy = 0; iy < nx[1]; iy++) {
    double yi = (double)iy * size[1];
    for (int iz = 0; iz < nx[2]; iz++) {
      double zi = (double)iz * size[2];
      double xi = FindX(section, yi, zi);

      fprintf(fp, "%e %e %e\n", xi, yi, zi);
    }
  }

  // セルを構成する点番号の書き出し
  fprintf(fp, "CELLS %d %d\n", nCells, 5 * nCells);
  for (int iy = 0; iy < nx[1] - 1; iy++) {
    for (int iz = 0; iz < nx[2] - 1; iz++) {
      int nodeID0 = nx[2] * (iy + 0) + iz + 0;
      int nodeID1 = nx[2] * (iy + 1) + iz + 0;
      int nodeID2 = nx[2] * (iy + 1) + iz + 1;
      int nodeID3 = nx[2] * (iy + 0) + iz + 1;

      fprintf(fp, "4 %d %d %d %d\n", nodeID0, nodeID1, nodeID2, nodeID3);
    }
  }

  // セルタイプの書き出し（9: Rectangular）
  fprintf(fp, "CELL_TYPES %d\n", nCells);
  for (int i = 0; i < nCells; i++) {
    fprintf(fp, "9\n");
  }

  // 応力の書き出し
  fprintf(fp, "POINT_DATA %d\n", nPoints);
  fprintf(fp, "SCALARS Stress float\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (int iy = 0; iy < nx[1]; iy++) {
    double yi = (double)iy * size[1];
    for (int iz = 0; iz < nx[2]; iz++) {
      double zi = (double)iz * size[2];
      double xi = FindX(section, yi, zi);
      double x[3] = {xi, yi, zi};

      fprintf(fp, "%e\n", Stress(volume, x));
    }
  }
}

static void XZSection(VOLUME_t *volume, FILE *fp) {
  SECTION_t *section = &volume->section;
  GRID_t *grid = &volume->grid;
  int nx[3] = {volume->nx[0] * (grid->nx[0] - 1) + 1,
               volume->nx[1] * (grid->nx[1] - 1) + 1,
               volume->nx[2] * (grid->nx[2] - 1) + 1};
  double *size = grid->size;
  int nPoints = nx[0] * nx[2];
  int nCells = (nx[0] - 1) * (nx[2] - 1);

  // 点の座標の書き出し
  fprintf(fp, "POINTS %d float\n", nPoints);
  for (int ix = 0; ix < nx[0]; ix++) {
    double xi = (double)ix * size[0];
    for (int iz = 0; iz < nx[2]; iz++) {
      double zi = (double)iz * size[2];
      double yi = FindY(section, xi, zi);

      fprintf(fp, "%e %e %e\n", xi, yi, zi);
    }
  }

  // セルを構成する点番号の書き出し
  fprintf(fp, "CELLS %d %d\n", nCells, 5 * nCells);
  for (int ix = 0; ix < nx[0] - 1; ix++) {
    for (int iz = 0; iz < nx[2] - 1; iz++) {
      int nodeID0 = nx[2] * (ix + 0) + iz + 0;
      int nodeID1 = nx[2] * (ix + 1) + iz + 0;
      int nodeID2 = nx[2] * (ix + 1) + iz + 1;
      int nodeID3 = nx[2] * (ix + 0) + iz + 1;

      fprintf(fp, "4 %d %d %d %d\n", nodeID0, nodeID1, nodeID2, nodeID3);
    }
  }

  // セルタイプの書き出し（9: Rectangular）
  fprintf(fp, "CELL_TYPES %d\n", nCells);
  for (int i = 0; i < nCells; i++) {
    fprintf(fp, "9\n");
  }

  // 応力の書き出し
  fprintf(fp, "POINT_DATA %d\n", nPoints);
  fprintf(fp, "SCALARS Stress float\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (int ix = 0; ix < nx[0]; ix++) {
    double xi = (double)ix * size[0];
    for (int iz = 0; iz < nx[2]; iz++) {
      double zi = (double)iz * size[2];
      double yi = FindY(section, xi, zi);
      double x[3] = {xi, yi, zi};

      fprintf(fp, "%e\n", Stress(volume, x));
    }
  }
}

static int SelectAxis(SECTION_t *section) {
  int axis = 0;
  double nmax = fabs(section->n[0]);

  // 法線ベクトルの成分の絶対値の最大値から基準となる軸を見つける
  if (nmax < fabs(section->n[1])) {
    nmax = fabs(section->n[1]);
    axis = 1;
  }
  if (nmax < fabs(section->n[2])) {
    nmax = fabs(section->n[2]);
    axis = 2;
  }
  return axis;
}

void WriteVolume(VOLUME_t *volume, char *directory) {
  char file_name[256];
  FILE *fp;
  double x0 = 0.0;
  double x1 = volume->size[0];
  double y0 = 0.0;
  double y1 = volume->size[1];
  double z0 = 0.0;
  double z1 = volume->size[2];

  sprintf(file_name, "%s/volume.vtk", directory);
  fp = fopen(file_name, "w");

  VTKHeader(fp);

  fprintf(fp, "POINTS %d float\n", 8);
  fprintf(fp, "%e %e %e\n", x0, y0, z0);
  fprintf(fp, "%e %e %e\n", x1, y0, z0);
  fprintf(fp, "%e %e %e\n", x1, y1, z0);
  fprintf(fp, "%e %e %e\n", x0, y1, z0);
  fprintf(fp, "%e %e %e\n", x0, y0, z1);
  fprintf(fp, "%e %e %e\n", x1, y0, z1);
  fprintf(fp, "%e %e %e\n", x1, y1, z1);
  fprintf(fp, "%e %e %e\n", x0, y1, z1);

  fprintf(fp, "CELLS %d %d\n", 1, 9);
  fprintf(fp, "8 0 1 2 3 4 5 6 7\n");

  fprintf(fp, "CELL_TYPES 1\n");
  fprintf(fp, "12\n");

  fclose(fp);
}

int main(int argc, char **argv) {
  VOLUME_t volume;
  FILE *fp;
  char file_name[256];

  // 入力データの読み込み
  ReadVolume(&volume, argv[1]);
  ReadStress(&volume, argv[1]);

  // 出力ファイルの準備
  sprintf(file_name, "%s/section.vtk", argv[1]);
  fp = fopen(file_name, "w");

  // VTKファイルのヘッダー部分の書き込み
  VTKHeader(fp);
  // 断面の法線ベクトルの成分の大きさから基準となる平面を選択する
  switch (SelectAxis(&volume.section)) {
    case 0:
      YZSection(&volume, fp);
      break;
    case 1:
      XZSection(&volume, fp);
      break;
    default:
      XYSection(&volume, fp);
      break;
  }
  fclose(fp);

  // ボリュームのVTKも書き出しておく
  WriteVolume(&volume, argv[1]);

  return 0;
}