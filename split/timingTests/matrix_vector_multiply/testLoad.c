#include "splitLoad.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
  if (argc != 2) {printf("Usage : testLoad filename"); exit(-1);}
  Data *dat = splitLoad(argv[1]);
  if (dat == NULL) {printf("Could not load data file"); exit(-1);}
  printData(dat, 0);

  Vector *x = getVar(dat, "x3");
  if (x == NULL) {printf("Could not find vector %s\n", "x3"); exit(-1);}
  printf("Got vec\n");

  printVector(x, 1);
  freeData(dat);
  return 0;
}
