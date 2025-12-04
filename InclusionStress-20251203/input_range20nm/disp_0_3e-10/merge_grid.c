#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  int n = atoi(argv[1]);
  int nx[3];
  int n_subgrid, nx_subgrid[3];
  FILE **ifp = (FILE **)malloc(n * sizeof(FILE *));
  FILE *ofp;
  char file_name[256];

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <number_of_files> <working_directory>\n",
            argv[0]);
    return 2;
  }

  for (int i = 0; i < n; i++) {
    sprintf(file_name, "%s/grid_stress_%04d.inp", argv[2], i + 1);
    ifp[i] = fopen(file_name, "r");
    if (ifp[i] == NULL) {
      fprintf(stderr, "Cannot open %s file\n", file_name);
      return 1;
    }
  }
  sprintf(file_name, "%s/grid_stress.inp", argv[2]);
  ofp = fopen(file_name, "w");

  for (int i = 0; i < n; i++) {
    fread(nx, sizeof(int), 3, ifp[i]);
  }
  fwrite(nx, sizeof(int), 3, ofp);

  for (int iGrid = 0; iGrid < nx[0] * nx[1] * nx[2]; iGrid++) {
    double stress[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (int i = 0; i < n; i++) {
      double istress[9];

      fread(istress, sizeof(double), 9, ifp[i]);
      for (int j = 0; j < 9; j++) {
        stress[j] += istress[j];
      }
    }
    fwrite(stress, sizeof(double), 9, ofp);
  }

  for (int i = 0; i < n; i++) {
    fread(&n_subgrid, sizeof(int), 1, ifp[i]);
  }
  fwrite(&n_subgrid, sizeof(int), 1, ofp);

  if (n_subgrid != 1) {
    for (int i = 0; i < n; i++) {
      fread(nx_subgrid, sizeof(int), 3, ifp[i]);
    }
    fwrite(nx_subgrid, sizeof(int), 3, ofp);
  }

  for (int iGrid = 0; iGrid < n_subgrid; iGrid++) {
    int grid_id;

    for (int i = 0; i < n; i++) {
      fread(&grid_id, sizeof(int), 1, ifp[i]);
    }
    fwrite(&grid_id, sizeof(int), 1, ofp);

    for (int jGrid = 0; jGrid < nx_subgrid[0] * nx_subgrid[1] * nx_subgrid[2]; jGrid++) {
      double stress[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      for (int i = 0; i < n; i++) {
        double istress[9];

        fread(istress, sizeof(double), 9, ifp[i]);
        for (int j = 0; j < 9; j++) {
          stress[j] += istress[j];
        }
      }
      fwrite(stress, sizeof(double), 9, ofp);
    }
  }

  for (int i = 0; i < n; i++) {
    fclose(ifp[i]);
  }
  fclose(ofp);

  return 0;
}