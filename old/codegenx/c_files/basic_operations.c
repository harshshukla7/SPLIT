#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header file containing function prototypes
#include "../header_files/basic_operations.h"

/********************************************************
 *          add two vectors of arbitrary size           *
 ********************************************************/

void add_vectors(double *resulting_vector, double *temp_vector1, double *temp_vector2, int size)
{
    
    int i;
    for(i = 0;i < size; i++)
    {
        resulting_vector[i] = temp_vector1[i] + temp_vector2[i];
    }
}

/********************************************************
 * compute the absolute value and save it for further   *
 ********************************************************/

void compute_absolute_value(double *dest_array, double *source_array, int size)
{
    int i;
    for(i = 0;i < size; i++)
    {
        dest_array[i] = fabs(source_array[i]);
    }
}

/*********************************************************
 *  compare elements with zero, return positive elements *
 *********************************************************/

void compare_elements_greater_than(double *dest_array, double *array, int size)
{
    int i;
    for(i = 0; i < size; i++)
    {
        if(array[i] > 0)
            dest_array[i] = 1;
        else
            dest_array[i] = 0;
    }
}

/*********************************************************
 *            compare elements: return 0 or 1            *
 *********************************************************/

void compare_elements_less_than(double *dest_array, double *array, int size)
{
    int i;
    
    for(i = 0; i < size; i++)
    {
        if(array[i] < 0)
            dest_array[i] = 1;
        else
            dest_array[i] = 0;
    }
}

/********************************************************
 * compare function for qsort function in the standard  *
 * library/ header file stdlib.h
 ********************************************************/

int compare_function (const void * a, const void * b)
{
    return ( *(int *)a - *(int *)b );
}

/*********************************************************************
 *              compute 2-norm from an arbitrary vector            *
 *********************************************************************/

void compute_norm(double *final_value, double *vector, int size)
{
    int i;
    double sum = 0.0;
    for(i = 0;i < size; i++)
    {
        sum = sum + vector[i]*vector[i];
    }
    *final_value = sqrt(sum);
    //*final_value = sum;
}

/*********************************************************************
 *              compute 1-norm from an arbitrary vector            *
 *********************************************************************/

void compute_norm_one(double *final_value, double *vector, int size)
{
    int i;
    double sum = 0.0;
    for(i = 0;i < size; i++)
    {
        sum = sum + fabs(vector[i]);
    }
    *final_value = sum;
}

/********************************************************
 *     compute sum of elements coming from an array     *
 ********************************************************/

void compute_sum(double *final_sum, double *array, int size)
{
    double sum = 0.0;
    int i;
    
    for(i = 0;i < size;i++)
    {
        sum = sum + array[i];
    }
    
    *final_sum = sum;
}

/*********************************************************
 * copy values from an array starting from a spec. index *
 *********************************************************/

void copy_n_values(double *array1, double *array2, int from, int size, int external_counter)
{
    int i, counter;
    
    counter = external_counter;
    for(i = from; i < size; i++)
    {
        array1[counter++] = array2[i];
    }
}

/*********************************************************
 * copy values from an array starting from a spec. index *
 *********************************************************/

void copy_matrix(int rows, int cols, double source_matrix[rows][cols], double dest_matrix[rows][cols])
{
    int i, j;

    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            dest_matrix[i][j] = source_matrix[i][j];
        }
    }
}

/********************************************************
 *    save the previous values of the opt. variables    *
 ********************************************************/

void copy_values(double *array1, double *array2, int size)
{
    int i;
    for(i = 0;i < size;i++)
    {
        array1[i] = array2[i];
    }
}

/*********************************************************
 *       create array of one: projection onto a box      *
 *********************************************************/

void create_array_of_ones(double *array, int size)
{
    int i;
    
    for(i = 0; i < size; i++)
    {
        array[i] = 1;
    }
}

/*********************************************************
 *                 compute dot product                   *
 *********************************************************/

void dot_product(double *array1, double *array2, double *result, int size)
{
    int i;
    int sum = 0.0f;    

    for(i = 0;i < size; i++)
    {
        sum = sum + array1[i]*array2[i];
    }

    *result = sum;
}

/*******************************************************
 *     compute element wise product of two vectors     *
 *******************************************************/

void element_wise_vector_multiplication(double *dest_vector, double *source_vector1, double *source_vector2, int size)
{
    int i;
        
    for(i = 0; i < size; i++)
    {
        dest_vector[i] = source_vector1[i] * source_vector2[i];
    }
}

/********************************************************
 *  find element with minimal value in the given array  *
 ********************************************************/

void find_minimum_element(double *vector, double *lower_bound, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (vector[i] >= lower_bound[i])
        {
            vector[i] = lower_bound[i];
        }
    }
}

/********************************************************
 *         compare elements and find the maximum        *
 ********************************************************/

void find_maximal_element(double *vector, double *value1, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (vector[i] <= *value1)
        {
            vector[i] = *value1;
        }
    }
}

/*******************************************************
 *      divide a 1 by each element from a vector       *
 *******************************************************/

void invert_elements(double *vector, int size)
{
    int i;
        
    for(i = 0; i < size; i++)
    {
        vector[i] = 1/vector[i];
    }
}

/*********************************************************
 *                matrix vector product                  *
 *********************************************************/

void matrix_vector_product(int rows, int cols, double matrix[rows][cols], double *array, double *res_vec)
{
    int i, j;
    double sum = 0.0f;
    
    for(i = 0; i < rows; i++)
    {
        sum = 0.0;
        for(j = 0; j < cols; j++)
        {
            sum = sum + matrix[i][j] * array[j];
        }
        res_vec[i] = sum;
    }
}

/*********************************************************
 *            multiply a matrix by a scalar              *
 *********************************************************/

void matrix_scalar_product(int rows, int cols, double matrix[rows][cols], double scalar)
{
    int i, j;  

    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            matrix[i][j] = scalar * matrix[i][j];
        }
    }
}

/*****************************************************
 * print the vector to the screen                    *
 *****************************************************/

void print_the_vector(double *A, int size)
{
    int i;
    for(i = 0; i < size;i++)
    {
        printf("%f\n", A[i]);
    }
}

/*****************************************************
 * print the matrix to the screen                    *
 *****************************************************/

void print_the_matrix(int rows, int cols, double matrix[rows][cols])
{
    int i, j;
    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            printf("%f ", matrix[i][j]);
        }
    printf("\n");
    }
}

/*********************************************************
 *                scalar vector addition                 *
 *********************************************************/

void scalar_vector_addition(double *dest_array, double *array, double scalar, int size)
{
    int i;
    
    for(i = 0; i < size; i++)
    {
        dest_array[i] = array[i] + scalar;
    }
}

/*********************************************************
 *          substract scalar from a vector               *
 *********************************************************/

void scalar_vector_substraction(double *dest_array, double *array, double scalar, int size)
{
    int i;
    
    for(i = 0; i < size; i++)
    {
        dest_array[i] = array[i] - scalar;
    }
}

/******************************************************
 *          multiply a vector with a scalar            *
 *******************************************************/

void scalar_vector_product(double *vector, double scalar, int size)
{
    int i;
    for(i = 0;i < size; i++)
    {
        vector[i] = vector[i] * scalar;
    }
}

/********************************************************
 *    extract values from a vector of arbitrary size    *
 *    expand these comments                             *
 ********************************************************/

void select_values(double *resulting_array, double *array, int size1, int size2)
{
    int i, index;
    index = 0;
    for(i = size1; i < size2; i++)
    {
        resulting_array[index] = array[i];
        index++;
    }
}

/********************************************************
 *     subtract two vectors of arbitrary size           *
 ********************************************************/

void subtract_vectors(double *resulting_vector, double *temp_vector1, double *temp_vector2, int size)
{
    int i;
    for(i = 0;i < size; i++)
    {
        resulting_vector[i] = temp_vector1[i] - temp_vector2[i];
    }
}

/********************************************************
 *                vector matrix product                 *
 ********************************************************/

void vector_matrix_product(int rows, int cols, double matrix[rows][cols], double *array, double *res_vec)
{
    int i, j;  
    double sum = 0.0f;

    for(i = 0; i < cols; i++)
    {
        sum = 0.0;
        for(j = 0; j < rows; j++)
        {
            sum = sum + array[j] * matrix[j][i];
        }
        res_vec[i] = sum;
    }
}

/********************************************************
 *   divide each element of the vector with a scalar    *
 ********************************************************/

void vector_scalar_division(double *vector, double scalar, int size)
{
    int i;
    for(i = 0;i < size; i++)
        vector[i] = vector[i]/scalar;
}