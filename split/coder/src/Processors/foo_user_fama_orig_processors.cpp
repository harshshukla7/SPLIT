/* 
* SPLIT code generation
*/


#include "foo_data.h"
#include "user_ama.h"

void foo_user(  data_t_state0_in state0_in_int[STATE0_IN_LENGTH],
				data_t_primal0_in primal0_in_int[PRIMAL0_IN_LENGTH],
				data_t_dual0_in dual0_in_int[DUAL0_IN_LENGTH],
				data_t_tol_iterates_in tol_iterates_in_int[TOL_ITERATES_IN_LENGTH],
				data_t_primal_out primal_out_int[PRIMAL_OUT_LENGTH],
				data_t_dual_out dual_out_int[DUAL_OUT_LENGTH],
				data_t_aux_primal_out aux_primal_out_int[AUX_PRIMAL_OUT_LENGTH],
				data_t_iterates_out iterates_out_int[ITERATES_OUT_LENGTH])
{

	Sol sol;
    int number_of_solves = 1;
    Opt opt = {tol_iterates_in_int[0],tol_iterates_in_int[1], tol_iterates_in_int[2], tol_iterates_in_int[3]};
    
    initialize();
    
    solve(&sol, par_ex, &opt);
    
    //aux_primal_out_int[0] = state0_in_int[0];
    
    copy_vector(primal_out_int, sol.primal, PRIMAL_OUT_LENGTH);
    copy_vector(dual_out_int, sol.dual, DUAL_OUT_LENGTH);
    copy_vector(aux_primal_out_int, sol.aux_prim, AUX_PRIMAL_OUT_LENGTH);
    
    iterates_out_int[2] = sol.itr;
    iterates_out_int[1] = sol.rDual;
    iterates_out_int[0] = sol.rPrimal;
    
}
