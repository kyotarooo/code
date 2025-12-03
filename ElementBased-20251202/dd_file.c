#include <stdio.h>

#include "dd_struct.h"

void DD_OpenLogFiles(DD_t *dd, char *directory) {
  char file_name[256];
  FILE *fp;
  FILE_t *file = &dd->file;

  // Direct interaction log
  // Check if the restart data file exsits
  sprintf(file_name, "%s/restart.dat", directory);
  fp = fopen(file_name, "r");
  sprintf (file_name, "%s/direct_interaction.log", directory);
  if (fp == NULL) {
    file->direct_interaction = fopen(file_name, "w");
  } else {
    file->direct_interaction = fopen(file_name, "a");
    fclose(fp);
  }
}