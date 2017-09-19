/* 
* icl::protoip
* Authors: asuardi <https://github.com/asuardi>, bulatkhusainov <https://github.com/bulatkhusainov>
* Date: November - 2014
*/


#include "soc_user.h"
#include "foo_function_wrapped.h"


void soc_user(float soc_state0_in[SOC_STATE0_IN_VECTOR_LENGTH],float soc_primal0_in[SOC_PRIMAL0_IN_VECTOR_LENGTH],float soc_dual0_in[SOC_DUAL0_IN_VECTOR_LENGTH],float soc_tol_iterates_in[SOC_TOL_ITERATES_IN_VECTOR_LENGTH],float soc_primal_out[SOC_PRIMAL_OUT_VECTOR_LENGTH],float soc_dual_out[SOC_DUAL_OUT_VECTOR_LENGTH],float soc_aux_primal_out[SOC_AUX_PRIMAL_OUT_VECTOR_LENGTH],float soc_iterates_out[SOC_ITERATES_OUT_VECTOR_LENGTH])
{
	// declare input and output interfaces arrays
	float *pl_vec_in_in,*pl_mat_in_in;
	float *pl_vec_out_out;


	// send data to DDR
	send_pl_vec_in_in(pl_vec_in_in);
	send_pl_mat_in_in(pl_mat_in_in);


	// call hardware accelerator assuming all interfaces are involoved
	start_foo(1,1,1);


	// wait for IP to finish
	while(!(finished_foo())){;}


	// read data from DDR
	receive_pl_vec_out_out(pl_vec_out_out);


}
