#include <Accelerate/Accelerate.h>
#include <stdio.h>

float a[800][800];

float x[800];// = {1.0f, 2.0f, 4.0f, 8.0f};  // the vector to be multiplied
float y[800];
 // = {0.f, 0.f, 0.f, 0.f,       // the result vector
 //    0.f, 0.f, 0.f, 0.f};


int main() {
    int i;

    for(int i=0; i<1000; i++) {
      cblas_sgemv(CblasRowMajor, CblasNoTrans, 800, 800, 1.0f, (float*)a, 800, x, 1, 1.0f, y, 1);
    }

    // for (i = 0; i < 8; i++) {
    //     printf("%.4f\n", y[i]);
    // }

    printf("Done!\n");

    return 0;
}