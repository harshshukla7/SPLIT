/* Load data file for SPLIT optimizer */

#ifndef __SPLITLOAD_H__
#define __SPLITLOAD_H__


// Data variable
#define CHARS_IN_NAME 100
struct Vector {
  unsigned int type;        // 1 => int, 0 => double
  unsigned int len;         // Number of elements
  char name[CHARS_IN_NAME]; // Name of the variable
  void *data;               // Vector of sizeof(type)*len
};
typedef struct Vector Vector;

struct Data {
  char info[255];          // Text string describing the file
  unsigned int numVectors; // Number of elements in vec
  Vector *vec;             // Holds the loaded variables
};
typedef struct Data Data;

// Load the data from filename
// Returns the variables from the file, or NULL if error
Data* splitLoad(char *filename);

// Gets the variable 'var' from the data set, or NULL if not found
Vector* getVar(Data *dat, char *var);

// Copy the variable 'var' to the pre-allocated array x.
// If var doesn't exist, it does nothing
void copyVar(void *x, Data *dat, char *var);

// Print out the loaded data
// verbosity = 0 : basic info
// verbosity > 1 : print variable values too
void printData(Data *dat, int verbosity);
void printVector(Vector *vec, int verbosity);

// Print out an unstructured vector
void printDoubleVec(double *vec, char *name, int len); 
void printFloatVec(float *vec, char *name, int len); 
void printIntVec(int *vec, char *name, int len); 

// Free memory for dat
void freeData(Data *dat);



#endif
