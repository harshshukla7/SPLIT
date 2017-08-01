/* Load data file for SPLIT optimizer */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "user_splitLoad.h"

#define DEBUG 0

#define TYPE_INT    0
#define TYPE_DOUBLE 2
#define TYPE_FLOAT  1

void* error(char *msg) {
  // TODO - this should also cleanup any allocated memory before returning
  printf("<ERROR> %s\n", msg);
  return NULL;
}

// Load the data from filename
// Returns the variables from the file
Data* splitLoad(char *filename) {
  FILE *f = fopen(filename, "r");
  if( f == NULL ) error("Could not open data file");

  Data *hdr = (Data*)malloc(sizeof(Data));
  size_t num_read;

  if (1   != fread(&(hdr->numVectors), sizeof(unsigned int), 1, f)) return error("Could not read data file when trying to get number of vectors");
  if (255 != fread(hdr->info,          sizeof(char), 255, f))       return error("Could not read data file when reading info");

  hdr->vec = (Vector*) malloc(sizeof(Vector) * hdr->numVectors);

  for (int i=0; i<hdr->numVectors; i++) {
    if (DEBUG) printf("Reading vector number %i", i);

    // Read the header
    // We do it this way because the compiler will change the size of structures by adding padding bytes
    if (1             != fread(&(hdr->vec[i].type), sizeof(unsigned int), 1, f))       return error("Could not read data file 1");
    if (1             != fread(&(hdr->vec[i].len),  sizeof(unsigned int), 1, f))       return error("Could not read data file 2");
    if (CHARS_IN_NAME != fread(hdr->vec[i].name,    sizeof(char),   CHARS_IN_NAME, f)) return error("Could not read data file 3");

    if (DEBUG) printf(" : name %s, len %i, type %i\n", hdr->vec[i].name, hdr->vec[i].len, hdr->vec[i].type);

    // Allocate memory
    size_t sz;
    switch(hdr->vec[i].type) {
      case TYPE_DOUBLE: sz = sizeof(double); break;
      case TYPE_INT:    sz = sizeof(int);    break;
      case TYPE_FLOAT:  sz = sizeof(float);  break;
      default: return error("Unknown variable type");
    }
    hdr->vec[i].data = malloc(hdr->vec[i].len * sz);
    if (hdr->vec[i].data == NULL) return error("Could not read variable header");

    if (hdr->vec[i].len != fread(hdr->vec[i].data,    sz, hdr->vec[i].len, f)) return error("Could not read data file 4");
  }

  fclose(f);

  return hdr;
}

// Print out the loaded data
// verbosity = 0 : basic info
// verbosity > 1 : print variable values too
void printData(Data *dat, int verbosity) {
  printf("%s\n", dat->info);
  printf("%i variables in data set\n", dat->numVectors);

  for (int i=0; i<dat->numVectors; i++) {
    printf("Var #%-3i : ", i);
    printVector(&(dat->vec[i]), verbosity);
  }
}

void printVector(Vector *vec, int verbosity) {
  printf("%-20s (len = %-4i, type = ", vec->name, vec->len);
  switch (vec->type) 
  {
    case TYPE_INT:    printf("int");    break;
    case TYPE_FLOAT:  printf("float");  break;
    case TYPE_DOUBLE: printf("double"); break;
  }
  printf(")\n");

  if (verbosity > 0) {
    printf("\t[");
    switch (vec->type) 
    {
      case TYPE_INT:    
      for (int j=0; j<vec->len; j++) printf("%i, ", ((int*) vec->data)[j]  );
        break;
      case TYPE_FLOAT:
      for (int j=0; j<vec->len; j++) printf("%.4f, ", ((float*) vec->data)[j]  );
        break;
      case TYPE_DOUBLE:
      for (int j=0; j<vec->len; j++) printf("%.4f, ", ((double*) vec->data)[j]  );
        break;
    }
    printf("]\n");
  }
}

// Print out an unstructured vector
void printDoubleVec(double *vec, char *name, int len) {
  printf("%-20s (len = %-4i, type = double)\n", name, len);
  printf("\t[");
  for (int j=0; j<len; j++) printf("%.4f, ", vec[j]);
    printf("]\n");
}
void printFloatVec(float *vec, char *name, int len) {
  printf("%-20s (len = %-4i, type = float)\n", name, len);
  printf("\t[");
    for (int j=0; j<len; j++) printf("%.4f, ", vec[j]);
      printf("]\n");
}
void printIntVec(int *vec, char *name, int len) {
  printf("%-20s (len = %-4i, type = int)\n", name, len);
  printf("\t[");
    for (int j=0; j<len; j++) printf("%i, ", vec[j]);
      printf("]\n");
}

// Free memory for dat
void freeData(Data *dat) {
  for (int i=0; i<dat->numVectors; i++) free(dat->vec[i].data);
    free(dat->vec);
  free(dat);
}

// Gets the variable 'var' from the data set, or NULL if error
Vector* getVar(Data *dat, char *var) {
  for (int i=0; i<dat->numVectors; i++) {
    if (strcmp(var, dat->vec[i].name) == 0) {
      return &(dat->vec[i]);
    }
  }
  return NULL;
}

// Copy the variable 'var' to the pre-allocated array x.
// If var doesn't exist, it does nothing
void copyVar(void *x, Data *dat, char *var) {
  int sz = 0;
  for (int i=0; i<dat->numVectors; i++) {
    if (strcmp(var, dat->vec[i].name) == 0) {
      switch (dat->vec[i].type) {
        case TYPE_INT:    sz = sizeof(int); break;
        case TYPE_FLOAT:  sz = sizeof(float); break;
        case TYPE_DOUBLE: sz = sizeof(double); break;
      }
      memcpy(x, dat->vec[i].data, sz*(dat->vec[i].len));
    }
  }
}
