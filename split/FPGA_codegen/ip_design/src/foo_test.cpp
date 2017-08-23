/* 
* icl::protoip
* Author: asuardi <https://github.com/asuardi>
* Date: November - 2014
*/


#include "foo_data.h"



void foo	(	
				uint32_t byte_pl_vec_in_in_offset,
				uint32_t byte_pl_mat_in_in_offset,
				uint32_t byte_pl_vec_out_out_offset,
				volatile data_t_memory *memory_inout);


using namespace std;
#define BUF_SIZE 64

//Input and Output vectors base addresses in the virtual memory
#define pl_vec_in_IN_DEFINED_MEM_ADDRESS 0
#define pl_mat_in_IN_DEFINED_MEM_ADDRESS (PL_VEC_IN_IN_LENGTH)*4
#define pl_vec_out_OUT_DEFINED_MEM_ADDRESS (PL_VEC_IN_IN_LENGTH+PL_MAT_IN_IN_LENGTH)*4


int main()
{

	char filename[BUF_SIZE]={0};

    int max_iter;

	uint32_t byte_pl_vec_in_in_offset;
	uint32_t byte_pl_mat_in_in_offset;
	uint32_t byte_pl_vec_out_out_offset;

	int32_t tmp_value;

	//assign the input/output vectors base address in the DDR memory
	byte_pl_vec_in_in_offset=pl_vec_in_IN_DEFINED_MEM_ADDRESS;
	byte_pl_mat_in_in_offset=pl_mat_in_IN_DEFINED_MEM_ADDRESS;
	byte_pl_vec_out_out_offset=pl_vec_out_OUT_DEFINED_MEM_ADDRESS;

	//allocate a memory named address of uint32_t or float words. Number of words is 1024 * (number of inputs and outputs vectors)
	data_t_memory *memory_inout;
	memory_inout = (data_t_memory *)malloc((PL_VEC_IN_IN_LENGTH+PL_MAT_IN_IN_LENGTH+PL_VEC_OUT_OUT_LENGTH)*4); //malloc size should be sum of input and output vector lengths * 4 Byte

	FILE *stimfile;
	FILE * pFile;
	int count_data;


	float *pl_vec_in_in;
	pl_vec_in_in = (float *)malloc(PL_VEC_IN_IN_LENGTH*sizeof (float));
	float *pl_mat_in_in;
	pl_mat_in_in = (float *)malloc(PL_MAT_IN_IN_LENGTH*sizeof (float));
	float *pl_vec_out_out;
	pl_vec_out_out = (float *)malloc(PL_VEC_OUT_OUT_LENGTH*sizeof (float));


	////////////////////////////////////////
	//read pl_vec_in_in vector

	// Open stimulus pl_vec_in_in.dat file for reading
	sprintf(filename,"pl_vec_in_in.dat");
	stimfile = fopen(filename, "r");

	// read data from file
	ifstream input1(filename);
	vector<float> myValues1;

	count_data=0;

	for (float f; input1 >> f; )
	{
		myValues1.push_back(f);
		count_data++;
	}

	//fill in input vector
	for (int i = 0; i<count_data; i++)
	{
		if  (i < PL_VEC_IN_IN_LENGTH) {
			pl_vec_in_in[i]=(float)myValues1[i];

			#if FLOAT_FIX_PL_VEC_IN_IN == 1
				tmp_value=(int32_t)(pl_vec_in_in[i]*(float)pow(2,(PL_VEC_IN_IN_FRACTIONLENGTH)));
				memory_inout[i+byte_pl_vec_in_in_offset/4] = *(uint32_t*)&tmp_value;
			#elif FLOAT_FIX_PL_VEC_IN_IN == 0
				memory_inout[i+byte_pl_vec_in_in_offset/4] = (float)pl_vec_in_in[i];
			#endif
		}

	}


	////////////////////////////////////////
	//read pl_mat_in_in vector

	// Open stimulus pl_mat_in_in.dat file for reading
	sprintf(filename,"pl_mat_in_in.dat");
	stimfile = fopen(filename, "r");

	// read data from file
	ifstream input2(filename);
	vector<float> myValues2;

	count_data=0;

	for (float f; input2 >> f; )
	{
		myValues2.push_back(f);
		count_data++;
	}

	//fill in input vector
	for (int i = 0; i<count_data; i++)
	{
		if  (i < PL_MAT_IN_IN_LENGTH) {
			pl_mat_in_in[i]=(float)myValues2[i];

			#if FLOAT_FIX_PL_MAT_IN_IN == 1
				tmp_value=(int32_t)(pl_mat_in_in[i]*(float)pow(2,(PL_MAT_IN_IN_FRACTIONLENGTH)));
				memory_inout[i+byte_pl_mat_in_in_offset/4] = *(uint32_t*)&tmp_value;
			#elif FLOAT_FIX_PL_MAT_IN_IN == 0
				memory_inout[i+byte_pl_mat_in_in_offset/4] = (float)pl_mat_in_in[i];
			#endif
		}

	}


	/////////////////////////////////////
	// foo c-simulation
	
	foo(	
				byte_pl_vec_in_in_offset,
				byte_pl_mat_in_in_offset,
				byte_pl_vec_out_out_offset,
				memory_inout);
	
	
	/////////////////////////////////////
	// read computed pl_vec_out_out and store it as pl_vec_out_out.dat
	pFile = fopen ("pl_vec_out_out.dat","w+");

	for (int i = 0; i < PL_VEC_OUT_OUT_LENGTH; i++)
	{

		#if FLOAT_FIX_PL_VEC_OUT_OUT == 1
			tmp_value=*(int32_t*)&memory_inout[i+byte_pl_vec_out_out_offset/4];
			pl_vec_out_out[i]=((float)tmp_value)/(float)pow(2,(PL_VEC_OUT_OUT_FRACTIONLENGTH));
		#elif FLOAT_FIX_PL_VEC_OUT_OUT == 0
			pl_vec_out_out[i]=(float)memory_inout[i+byte_pl_vec_out_out_offset/4];
		#endif
		
		fprintf(pFile,"%f \n ",pl_vec_out_out[i]);

	}
	fprintf(pFile,"\n");
	fclose (pFile);
		

	return 0;
}
