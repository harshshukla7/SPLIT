#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/basic_operations.h"
#include "../header_files/compute_array_xbar.h"

void compute_array_xbar(int x_len, double arrayX[x_len], double arrayXbar[x_len],
                        double arrayPrevX[x_len], double theta)

// function body                       
{
     // x_bar = x + theta*(x - x_prev)
    
     double temp_array[x_len];
     double temp_array1[x_len];
     
     // x - x_prev
     subtract_vectors(temp_array, arrayX, arrayPrevX, x_len);
     
     // theta*(x - x_prev)
     scalar_vector_product(temp_array, theta, x_len);
     
     // x + theta*(x - x_prev)
     add_vectors(temp_array1, arrayX, temp_array, x_len);
     copy_values(arrayXbar, temp_array1, x_len);
}