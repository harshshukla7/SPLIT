#include <stdio.h>

#include "splitLoad.h"
#include "splitData.h"

int main(void) {
  // Load the data manually, so that we can print it out
  Data* dat = splitLoad("splitData.dat");
  // Print out the variables
  printData(dat, 1);
  freeData(dat);  


  // Load the data the 'normal way'
  // The data will be in the variables defined in splitData.h
  loadData();
}