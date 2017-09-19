/* 
* icl::protoip
* Author: asuardi <https://github.com/asuardi>
* Date: November - 2014
*/


#include "foo_data.h"


void foo_user(  data_t_pl_vec_in_in pl_vec_in_in_int[PL_VEC_IN_IN_LENGTH],
				data_t_pl_mat_in_in pl_mat_in_in_int[PL_MAT_IN_IN_LENGTH],
				data_t_pl_vec_out_out pl_vec_out_out_int[PL_VEC_OUT_OUT_LENGTH])
{

	///////////////////////////////////////
	//ADD USER algorithm here below:
	//(this is an example)
	alg_0 : for(int i = 0; i <PL_VEC_OUT_OUT_LENGTH; i++)
	{
		pl_vec_out_out_int[i]=0;
		loop_0_pl_vec_in : for(int i_pl_vec_in = 0; i_pl_vec_in <PL_VEC_IN_IN_LENGTH; i_pl_vec_in++)
		{
			pl_vec_out_out_int[i]=pl_vec_out_out_int[i] + (data_t_pl_vec_out_out)pl_vec_in_in_int[i_pl_vec_in];
		}
		loop_0_pl_mat_in : for(int i_pl_mat_in = 0; i_pl_mat_in <PL_MAT_IN_IN_LENGTH; i_pl_mat_in++)
		{
			pl_vec_out_out_int[i]=pl_vec_out_out_int[i] + (data_t_pl_vec_out_out)pl_mat_in_in_int[i_pl_mat_in];
		}
	}

}
