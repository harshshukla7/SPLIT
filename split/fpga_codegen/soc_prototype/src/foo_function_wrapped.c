/* 
* icl::protoip
* Author: asuardi <https://github.com/asuardi>
* Date: November - 2014
*/


#include "foo_function_wrapped.h"

typedef uint32_t           Xint32;

XFoo xcore;

//functions for sending data from PS to DDR
void send_pl_vec_in_in(float* pl_vec_in_in)
{
	Xint32 *pl_vec_in_in_ptr_ddr = (Xint32 *)pl_vec_in_IN_DEFINED_MEM_ADDRESS;
	int32_t inputvec_fix[PL_VEC_IN_IN_VECTOR_LENGTH];
	int i;

	//write pl_vec_in_in to DDR
	if (FLOAT_FIX_PL_VEC_IN_IN == 1)
	{
		for(i = 0; i < PL_VEC_IN_IN_VECTOR_LENGTH; i++) // convert floating point to fixed
		{
			inputvec_fix[i] = (int32_t)(pl_vec_in_in[i]*pow(2, PL_VEC_IN_IN_FRACTIONLENGTH));
		}
		memcpy(pl_vec_in_in_ptr_ddr, inputvec_fix, PL_VEC_IN_IN_VECTOR_LENGTH*4);
	}
	else { //floating point
		memcpy(pl_vec_in_in_ptr_ddr, pl_vec_in_in, PL_VEC_IN_IN_VECTOR_LENGTH*4);
	}
}
void send_pl_mat_in_in(float* pl_mat_in_in)
{
	Xint32 *pl_mat_in_in_ptr_ddr = (Xint32 *)pl_mat_in_IN_DEFINED_MEM_ADDRESS;
	int32_t inputvec_fix[PL_MAT_IN_IN_VECTOR_LENGTH];
	int i;

	//write pl_mat_in_in to DDR
	if (FLOAT_FIX_PL_MAT_IN_IN == 1)
	{
		for(i = 0; i < PL_MAT_IN_IN_VECTOR_LENGTH; i++) // convert floating point to fixed
		{
			inputvec_fix[i] = (int32_t)(pl_mat_in_in[i]*pow(2, PL_MAT_IN_IN_FRACTIONLENGTH));
		}
		memcpy(pl_mat_in_in_ptr_ddr, inputvec_fix, PL_MAT_IN_IN_VECTOR_LENGTH*4);
	}
	else { //floating point
		memcpy(pl_mat_in_in_ptr_ddr, pl_mat_in_in, PL_MAT_IN_IN_VECTOR_LENGTH*4);
	}
}

//function for calling foo_user IP
void start_foo(uint32_t pl_vec_in_in_required,uint32_t pl_mat_in_in_required,uint32_t pl_vec_out_out_required)
{
	xcore.Bus_a_BaseAddress = 0x43c00000;
	xcore.IsReady = XIL_COMPONENT_IS_READY;

		if(pl_vec_in_in_required)
		{
			XFoo_Set_byte_pl_vec_in_in_offset(&xcore,pl_vec_in_IN_DEFINED_MEM_ADDRESS);
		}
		else
		{
			XFoo_Set_byte_pl_vec_in_in_offset(&xcore,(1<<31));
		}
		if(pl_mat_in_in_required)
		{
			XFoo_Set_byte_pl_mat_in_in_offset(&xcore,pl_mat_in_IN_DEFINED_MEM_ADDRESS);
		}
		else
		{
			XFoo_Set_byte_pl_mat_in_in_offset(&xcore,(1<<31));
		}
		if(pl_vec_out_out_required)
		{
			XFoo_Set_byte_pl_vec_out_out_offset(&xcore,pl_vec_out_OUT_DEFINED_MEM_ADDRESS);
		}
		else
		{
			XFoo_Set_byte_pl_vec_out_out_offset(&xcore,(1<<31));
		}
		XFoo_Start(&xcore);
}

//function for checking foo_user IP
uint32_t finished_foo(void)
{
	return XFoo_IsIdle(&xcore);
}
//functions for receiving data from DDR to PS
void receive_pl_vec_out_out(float* pl_vec_out_out)
{
	Xint32 *pl_vec_out_out_ptr_ddr = (Xint32 *)pl_vec_out_OUT_DEFINED_MEM_ADDRESS;
	int32_t outputvec_fix[PL_VEC_OUT_OUT_VECTOR_LENGTH];
	int i;

	//read pl_vec_out_out from DDR
	if (FLOAT_FIX_PL_VEC_OUT_OUT == 1) { //fixed point
		memcpy(outputvec_fix,pl_vec_out_out_ptr_ddr,  PL_VEC_OUT_OUT_VECTOR_LENGTH*4);
		for(i = 0; i < PL_VEC_OUT_OUT_VECTOR_LENGTH; i++)
		{
			pl_vec_out_out[i] = ((float)outputvec_fix[i]/pow(2, PL_VEC_OUT_OUT_FRACTIONLENGTH));
		}
	} else { //floating point
		memcpy(pl_vec_out_out,pl_vec_out_out_ptr_ddr,  PL_VEC_OUT_OUT_VECTOR_LENGTH*4);
	}
}
//--------------------------------------------------------------------
