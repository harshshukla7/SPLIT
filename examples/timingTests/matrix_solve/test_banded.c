#include <stdio.h>

#include "splitTimer.h"
#include "splitLoad.h"
#include <Accelerate/Accelerate.h>

#include "splitData.h"

int main() {
  FILE *f = fopen("data.dat", "a+");

  loadData();

  double x[N];
  test_banded(x,b);

  printVec(x, "x", N); 
}
