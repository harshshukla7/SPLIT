#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../header_files/projections.h"
#include "../header_files/basic_operations.h"

// compute projecion onto a box
void proj_box(double *x, double *z, double *lb, double *ub, int size)
{
    
    double temp_vector[size];
      
    // copy values for finding the maximal element
    copy_values(temp_vector, z, size);
    
    // compute maximal value
    find_maximal_element(temp_vector, lb, size);
    
    // compute maximal value
    find_minimum_element(temp_vector, ub, size);
        
    //compute value for finding the minimal value
    copy_values(x, temp_vector, size);

}

// compute projecion onto a nonnegative orthant
void proj_positive(double *x, double *z, int size, double weight, char *norm_type, double constant)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (z[i] <= 0)
            x[i] = 0;
        else
            x[i] = z[i];
    }   
}

// compute projecion onto a nonpositive orthant
void proj_negative(double *x, double *z, int size, double weight, char *norm_type, double constant)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (z[i] >= 0)
            x[i] = 0;
        else
            x[i] = z[i];
    }   
}

// compute projection - normBall
void proj_normBall(double *x, double *z, int size, double weight, char *norm_type, double constant)
{
     // variable declarations	
     int i, index_after_break;
     double nz_norm, final_sum, v_array_chunk_sum;
     double lam; // auxiliary variable needed in this particular algorithm
     double v_array[size];
     double v_array_temp[size];
     double v_array_final[size];
     double v_array_chunk[size]; // this is a bit memory wasting - reimplement it if it is possible
     double temp_array1[size], temp_array2[size], temp_array3[size], temp_array4[size];
     double aux_array1[size], aux_array2[size];
     
     // fix it later
     int const_value = 0;

     // lower and upper bound - infinity norm
     double pos_ones[size];
     double neg_ones[size];

     if (!strcmp(norm_type, "1"))
     {    
   	  // compute norm of vector z
    	  compute_norm_one(&nz_norm, z, size);
    
    	  /* compute projection */
          // Non-unitary stuff is not implemented yet, hence the value is fixed to be equal to 1
    	  if(nz_norm <= 1) 

	       lam = 0.0;

      	  else
          {
         	// lam solves the equation sum max(|x|-lam, 0) == 1
         	compute_absolute_value(v_array, z, size);
         	qsort(v_array, size, sizeof(double), compare_function);
        
         	for(i = 1;i <= size;i++)
         	{
            	  scalar_vector_substraction(v_array_temp, v_array, v_array[i-1], size);
                  proj_positive(v_array_final, v_array_temp, size, weight, norm_type, constant);
                  compute_sum(&final_sum, v_array_temp, size);
                  if(final_sum < 1 )
                  {
                    index_after_break = i - 1;
                    break;
                  } 
            }
         	// we know that lam lies between v(i-1) and v(i), and that i > 1
         	copy_n_values(v_array_chunk, v_array, index_after_break, size, const_value);
         	compute_sum(&v_array_chunk_sum, v_array_chunk, size - index_after_break);
         	v_array_chunk_sum = v_array_chunk_sum - 1;
         	lam = v_array_chunk_sum / (size - i + 1);
          }
 
         // z + lam
         scalar_vector_addition(temp_array1, z, lam, size);
        
         // z - lam
         scalar_vector_substraction(temp_array2, z, lam, size);
        
         // z + lam < 0
         compare_elements_less_than(temp_array3, temp_array1, size);
                
         // lam - z > 0
         compare_elements_greater_than(temp_array4, temp_array2, size);
        
         // compute projection, store the result in x
         element_wise_vector_multiplication(aux_array1, temp_array3, temp_array1, size);
         element_wise_vector_multiplication(aux_array2, temp_array4, temp_array2, size);
         add_vectors(x, aux_array1, aux_array2, size);

     } // end of the first if
    
     else if(!strcmp(norm_type, "2"))
     {
         
         // compute the square root rather here
         // It looks like I do not need it
         // constant = sqrt(constant);
         
    	  // compute norm of vector y
         compute_norm(&nz_norm, z, size);
       
         // evaluate projection - Je pense
         if(nz_norm <= constant)
         {
            for(i = 0; i < size; i++)
            {
              x[i] = z[i];
            }
         }
         else
         {
           for(i = 0;i < size;i++)
           {
              copy_values(x, z, size);
              scalar_vector_product(x, constant/nz_norm, size);
              // vector_scalar_division(x, nz_norm, size);
           }
         }
     } // end of elif

     else if(!strcmp(norm_type, "inf") || !strcmp(norm_type, "Inf"))
     {    
          // prepare lower bound and scale it properly
          create_array_of_ones(neg_ones, size);
          scalar_vector_product(neg_ones, -1*constant, size);
          
          // prepare upper bound and scale it properly
          create_array_of_ones(pos_ones, size);
          scalar_vector_product(pos_ones, constant, size);
                    
          // compute projection onto a box
	      proj_box(x, z, neg_ones, pos_ones, size); 
     
     }
     else
         fprintf(stderr, "Unknown norm");
}

// compute projection quad ball
void proj_quadball(double *x, double *z, int size, double weight, char *norm_type, double constant)
{
     double nz_norm;
     int i;
     
     // compute the square root rather here
     // I need it here: I think
     constant = sqrt(constant);
    
	 // compute norm of vector y
     compute_norm(&nz_norm, z, size);
        
     // evaluate projection - Je pense
     if(nz_norm <= constant)
     {
        for(i = 0; i < size; i++)
        {
           x[i] = z[i];
        }
     }
     else
     {
        for(i = 0;i < size;i++)
        {
           copy_values(x, z, size);
           scalar_vector_product(x, constant/nz_norm, size);
        }
     }
}

// compute projection soc
void proj_socp(double *x, double *z, int size, double weight, char *norm_type, double constant)
{
    double t;
    int i;
    double y[size - 1];
    double ny_norm, proj_const;
    
    // assign the first element from vector z into the variable t
    t = z[0];
        
    // assign the last n-1 element from vector z into the array y
    for(i = 0; i < size-1; i++)
    {
        y[i] = z[i+1];
    }
    
    // compute norm of vector y
    compute_norm(&ny_norm, y, size-1);
     
    // evaluate projection - Je pense
    if(ny_norm <= t)
    {
       for(i = 0; i < size; i++)
       {
           x[i] = z[i];
       }
    }
    else if(ny_norm <= -t)
    {
        for(i = 0;i < size;i++)
        {
            x[i] = 0;
        }
    }
    else
    {
        // scalar in the projection phase: Colin's code: (t + ny)/(2 * ny)
        proj_const = (t + ny_norm) / (2 * ny_norm);
        
        // prepare vector being multiplied by a scalar
        x[0] = ny_norm; // first element
           
        // remaining elements
        for(i = 1;i < size; i++)
        {
            x[i] = y[i-1];
        }
       
        scalar_vector_product(x, proj_const, size);
    }
}

// compute projection conjugate of soc
void proj_socp_conj(double *x, double *z, int size, double weight, char *norm_type, double constant)
{
    double t;
    int i;
    double y[size - 1];
    double ny_norm, proj_const;
    double z_copy[size], x_copy[size];
    
    // modify the input vector:
    // z = [-z(1); z(2:end)];
    
    copy_values(z_copy, z, size);
    z[0] = -1.0 * z_copy[0];
    copy_n_values(z, z_copy, 1, size, 1);
      
    // assign the first element from vector z into the variable t
    t = z[0];
        
    // assign the last n-1 element from vector z into the array y
    for(i = 0; i < size-1; i++)
    {
        y[i] = z[i+1];
    }
    
    // compute norm of vector y
    compute_norm(&ny_norm, y, size-1);
     
    // evaluate projection - Je pense
    if(ny_norm <= t)
    {
       for(i = 0; i < size; i++)
       {
           x[i] = z[i];
       }
    }
    else if(ny_norm <= -t)
    {
        for(i = 0;i < size;i++)
        {
            x[i] = 0;
        }
    }
    else
    {
        // scalar in the projection phase: Colin's code: (t + ny)/(2 * ny)
        proj_const = (t + ny_norm) / (2 * ny_norm);
        
        // prepare vector being multiplied by a scalar
        x[0] = ny_norm; // first element
           
        // remaining elements
        for(i = 1;i < size; i++)
        {
            x[i] = y[i-1];
        }
       
        scalar_vector_product(x, proj_const, size);
    }
    
    // modify the output vector
    // x = [-x(1,1); x(2:end,1)];
    copy_values(x_copy, x, size);
    x[0] = -1.0 * x_copy[0];
    copy_n_values(x, x_copy, 1, size, 1);
}

// compute prox operator associated with the corresponding norm
void prox_norm(double *x, double *z, int size, double weight, char *norm_type, double constant)
{
    char  pDual[4];      // dual to the corresponding norm
    double z_temp[size];  // vector z divided by scalar t
    double x_temp[size];  // 
    
    if (!strcmp(norm_type, "1"))
         strcpy(pDual, "Inf");
    else if(!strcmp(norm_type, "2"))
         strcpy(pDual, "2");
    else if(!strcmp(norm_type, "inf") || !strcmp(norm_type, "Inf"))
         strcpy(pDual, "1");
    else
         fprintf(stderr, "Unknown norm");
    
    //compute projection 
    copy_values(z_temp, z, size);
    vector_scalar_division(z_temp, weight, size);    
    proj_normBall(x_temp, z_temp, size, weight, pDual, 1);
    scalar_vector_product(x_temp, weight, size);
    subtract_vectors(x, z, x_temp, size);
}