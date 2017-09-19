/* 
* icl::protoip
* Authors: asuardi <https://github.com/asuardi>, bulatkhusainov <https://github.com/bulatkhusainov>
* Date: November - 2014
*/


#include "soc_user.h"
#include "foo_function_wrapped.h"
#include "user_ama.h"

void soc_user(float soc_state0_in[SOC_STATE0_IN_VECTOR_LENGTH],float soc_primal0_in[SOC_PRIMAL0_IN_VECTOR_LENGTH],float soc_dual0_in[SOC_DUAL0_IN_VECTOR_LENGTH],float soc_tol_iterates_in[SOC_TOL_ITERATES_IN_VECTOR_LENGTH],float soc_primal_out[SOC_PRIMAL_OUT_VECTOR_LENGTH],float soc_dual_out[SOC_DUAL_OUT_VECTOR_LENGTH],float soc_aux_primal_out[SOC_AUX_PRIMAL_OUT_VECTOR_LENGTH],float soc_iterates_out[SOC_ITERATES_OUT_VECTOR_LENGTH])
{
	// declare input and output interfaces arrays
			int i;

	    float balance;

	    Sol sol;
	    int number_of_solves = 1;
	    Opt opt = {soc_tol_iterates_in[0], soc_tol_iterates_in[1], soc_tol_iterates_in[2], soc_tol_iterates_in[3]};

	//copy_vector(par_ex, soc_state0_in,  SOC_STATE0_IN_VECTOR_LENGTH );
	initialize();
	solve(&sol, soc_state0_in, &opt);
	/*for( i=0;i<SOC_SOC_X_OUT_OUT_VECTOR_LENGTH;i++){

	//balance = soc_x_hat_in[i]*soc_x_hat_in[i];

	balance = (float) sol.primal[i];
	soc_soc_x_out_out[i] = balance;

	}*/

    copy_vector(soc_primal_out, sol.primal, SOC_PRIMAL_OUT_VECTOR_LENGTH);
    copy_vector(soc_dual_out, sol.dual, SOC_DUAL_OUT_VECTOR_LENGTH);
    copy_vector(soc_aux_primal_out, sol.aux_prim, SOC_AUX_PRIMAL_OUT_VECTOR_LENGTH);
    
    soc_iterates_out[2] = sol.itr;
    soc_iterates_out[1] = sol.rDual;
    soc_iterates_out[0] = sol.rPrimal;


}
