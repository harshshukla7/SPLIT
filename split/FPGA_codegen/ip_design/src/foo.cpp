/* 
* icl::protoip
* Author: asuardi <https://github.com/asuardi>
* Date: November - 2014
*/


#include "foo_data.h"


void foo_user(  data_t_pl_vec_in_in pl_vec_in_in_int[PL_VEC_IN_IN_LENGTH],
				data_t_pl_mat_in_in pl_mat_in_in_int[PL_MAT_IN_IN_LENGTH],
				data_t_pl_vec_out_out pl_vec_out_out_int[PL_VEC_OUT_OUT_LENGTH]);


void foo	(
				uint32_t byte_pl_vec_in_in_offset,
				uint32_t byte_pl_mat_in_in_offset,
				uint32_t byte_pl_vec_out_out_offset,
				volatile data_t_memory *memory_inout)
{

	#ifndef __SYNTHESIS__
	//Any system calls which manage memory allocation within the system, for example malloc(), alloc() and free(), must be removed from the design code prior to synthesis. 

	data_t_interface_pl_vec_in_in *pl_vec_in_in;
	pl_vec_in_in = (data_t_interface_pl_vec_in_in *)malloc(PL_VEC_IN_IN_LENGTH*sizeof(data_t_interface_pl_vec_in_in));
	data_t_interface_pl_mat_in_in *pl_mat_in_in;
	pl_mat_in_in = (data_t_interface_pl_mat_in_in *)malloc(PL_MAT_IN_IN_LENGTH*sizeof(data_t_interface_pl_mat_in_in));
	data_t_interface_pl_vec_out_out *pl_vec_out_out;
	pl_vec_out_out = (data_t_interface_pl_vec_out_out *)malloc(PL_VEC_OUT_OUT_LENGTH*sizeof(data_t_interface_pl_vec_out_out));

	data_t_pl_vec_in_in *pl_vec_in_in_int;
	pl_vec_in_in_int = (data_t_pl_vec_in_in *)malloc(PL_VEC_IN_IN_LENGTH*sizeof (data_t_pl_vec_in_in));
	data_t_pl_mat_in_in *pl_mat_in_in_int;
	pl_mat_in_in_int = (data_t_pl_mat_in_in *)malloc(PL_MAT_IN_IN_LENGTH*sizeof (data_t_pl_mat_in_in));
	data_t_pl_vec_out_out *pl_vec_out_out_int;
	pl_vec_out_out_int = (data_t_pl_vec_out_out *)malloc(PL_VEC_OUT_OUT_LENGTH*sizeof (data_t_pl_vec_out_out));

	#else
	//for synthesis

	data_t_interface_pl_vec_in_in  pl_vec_in_in[PL_VEC_IN_IN_LENGTH];
	data_t_interface_pl_mat_in_in  pl_mat_in_in[PL_MAT_IN_IN_LENGTH];
	data_t_interface_pl_vec_out_out  pl_vec_out_out[PL_VEC_OUT_OUT_LENGTH];

	static data_t_pl_vec_in_in  pl_vec_in_in_int[PL_VEC_IN_IN_LENGTH];
	static data_t_pl_mat_in_in  pl_mat_in_in_int[PL_MAT_IN_IN_LENGTH];
	data_t_pl_vec_out_out  pl_vec_out_out_int[PL_VEC_OUT_OUT_LENGTH];

	#endif

	#if FLOAT_FIX_PL_VEC_IN_IN == 1
	///////////////////////////////////////
	//load input vectors from memory (DDR)

	if(!(byte_pl_vec_in_in_offset & (1<<31)))
	{
		memcpy(pl_vec_in_in,(const data_t_memory*)(memory_inout+byte_pl_vec_in_in_offset/4),PL_VEC_IN_IN_LENGTH*sizeof(data_t_memory));

    	//Initialisation: cast to the precision used for the algorithm
		input_cast_loop_pl_vec_in:for (int i=0; i< PL_VEC_IN_IN_LENGTH; i++)
			pl_vec_in_in_int[i]=(data_t_pl_vec_in_in)pl_vec_in_in[i];

	}
	

	#elif FLOAT_FIX_PL_VEC_IN_IN == 0
	///////////////////////////////////////
	//load input vectors from memory (DDR)

	if(!(byte_pl_vec_in_in_offset & (1<<31)))
	{
		memcpy(pl_vec_in_in_int,(const data_t_memory*)(memory_inout+byte_pl_vec_in_in_offset/4),PL_VEC_IN_IN_LENGTH*sizeof(data_t_memory));
	}

	#endif


	#if FLOAT_FIX_PL_MAT_IN_IN == 1
	///////////////////////////////////////
	//load input vectors from memory (DDR)

	if(!(byte_pl_mat_in_in_offset & (1<<31)))
	{
		memcpy(pl_mat_in_in,(const data_t_memory*)(memory_inout+byte_pl_mat_in_in_offset/4),PL_MAT_IN_IN_LENGTH*sizeof(data_t_memory));

    	//Initialisation: cast to the precision used for the algorithm
		input_cast_loop_pl_mat_in:for (int i=0; i< PL_MAT_IN_IN_LENGTH; i++)
			pl_mat_in_in_int[i]=(data_t_pl_mat_in_in)pl_mat_in_in[i];

	}
	

	#elif FLOAT_FIX_PL_MAT_IN_IN == 0
	///////////////////////////////////////
	//load input vectors from memory (DDR)

	if(!(byte_pl_mat_in_in_offset & (1<<31)))
	{
		memcpy(pl_mat_in_in_int,(const data_t_memory*)(memory_inout+byte_pl_mat_in_in_offset/4),PL_MAT_IN_IN_LENGTH*sizeof(data_t_memory));
	}

	#endif



	///////////////////////////////////////
	//USER algorithm function (foo_user.cpp) call
	//Input vectors are:
	//pl_vec_in_in_int[PL_VEC_IN_IN_LENGTH] -> data type is data_t_pl_vec_in_in
	//pl_mat_in_in_int[PL_MAT_IN_IN_LENGTH] -> data type is data_t_pl_mat_in_in
	//Output vectors are:
	//pl_vec_out_out_int[PL_VEC_OUT_OUT_LENGTH] -> data type is data_t_pl_vec_out_out
	foo_user_top: foo_user(	pl_vec_in_in_int,
							pl_mat_in_in_int,
							pl_vec_out_out_int);


	#if FLOAT_FIX_PL_VEC_OUT_OUT == 1
	///////////////////////////////////////
	//store output vectors to memory (DDR)

	if(!(byte_pl_vec_out_out_offset & (1<<31)))
	{
		output_cast_loop_pl_vec_out: for(int i = 0; i <  PL_VEC_OUT_OUT_LENGTH; i++)
			pl_vec_out_out[i]=(data_t_interface_pl_vec_out_out)pl_vec_out_out_int[i];

		//write results vector y_out to DDR
		memcpy((data_t_memory *)(memory_inout+byte_pl_vec_out_out_offset/4),pl_vec_out_out,PL_VEC_OUT_OUT_LENGTH*sizeof(data_t_memory));

	}
	#elif FLOAT_FIX_PL_VEC_OUT_OUT == 0
	///////////////////////////////////////
	//write results vector y_out to DDR
	if(!(byte_pl_vec_out_out_offset & (1<<31)))
	{
		memcpy((data_t_memory *)(memory_inout+byte_pl_vec_out_out_offset/4),pl_vec_out_out_int,PL_VEC_OUT_OUT_LENGTH*sizeof(data_t_memory));
	}

	#endif




}
