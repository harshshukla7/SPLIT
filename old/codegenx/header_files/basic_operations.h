#ifndef _BASIC_OPERATIONS__H__
#define _BASIC_OPERATIONS__H__


	void add_vectors(double *resulting_vector, double *temp_vector1, double *temp_vector2, int size);
    
	void compute_absolute_value(double *dest_array, double *source_array, int size);
	void compare_elements_greater_than(double *dest_array, double *array, int size);
	void compare_elements_less_than(double *dest_array, double *array, int size);
	int  compare_function (const void * a, const void * b);
	void compute_norm(double *final_value, double *vector, int size);
	void compute_norm_one(double *final_value, double *vector, int size);
	void compute_sum(double *final_sum, double *array, int size);
	void copy_n_values(double *array1, double *array2, int from, int size, int external_counter);
    void copy_matrix(int rows, int cols, double source_matrix[rows][cols], double dest_matrix[rows][cols]);
    void copy_values(double *array1, double *array2, int size);
    void create_array_of_ones(double *array, int size);

    void dot_product(double *array1, double *array2, double *result, int size);

	void element_wise_vector_multiplication(double *dest_vector, double *source_vector1, double *source_vector2, int size);

	void find_maximal_element(double *vector, double *value1, int size);
	void find_minimum_element(double *vector, double *lower_bound, int size);
    
    void invert_elements(double *vector, int size);

    void matrix_vector_product(int rows, int cols, double matrix[rows][cols], double *array, double *res_vec);
    void matrix_scalar_product(int rows, int cols, double matrix[rows][cols], double scalar);
        
	void print_the_vector(double *A, int size);
    void print_the_matrix(int rows, int cols, double matrix[rows][cols]);

	void scalar_vector_addition(double *dest_array, double *array, double scalar, int size);
	void scalar_vector_substraction(double *dest_array, double *array, double scalar, int size);
	void scalar_vector_product(double *vector, double scalar, int size);
    void select_values(double *resulting_array, double *array, int size1, int size2);
    void subtract_vectors(double *resulting_vector, double *temp_vector1, double *temp_vector2, int size);

    void vector_matrix_product(int rows, int cols, double matrix[rows][cols], double *array, double *res_vec);	
	void vector_scalar_division(double *vector, double scalar, int size);

#endif