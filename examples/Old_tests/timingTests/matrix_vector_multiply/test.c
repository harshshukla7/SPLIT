#include <stdio.h>
#include <string.h>

float sum_vector(int n, float *v) { int i; float s = 0.0; for(i=0; i<n; i++) { s = s + v[i]; } return s; }

int main(void) {
  float v[10000];
  int   n = 10000;

  for(int i=0; i<100000; i++)
    sum_vector(n,v);
  printf("Hello world\n");

  return 0;
}

// void test1(int *List, int Length) {
//   int i = 0;
//   while(i < Length) {
//     List[i] = i*2;
//     i++;
//   }
// }

// int main(void) {
//   int List[40];
//   int Length = 40;

//   test1(List, Length);

//   printf("Hello world\n");

//   return 0;
// }
