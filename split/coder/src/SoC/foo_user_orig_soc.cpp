
/* 
* icl::protoip-SPLIT
*/


#include "foo_data.h"
#include "user_mv_mult.h" 

void foo_user(  data_t_pl_vec_in_in pl_vec_in_in_int[PL_VEC_IN_IN_LENGTH],
				data_t_pl_vec_out_out pl_vec_out_out_int[PL_VEC_OUT_OUT_LENGTH])
{


    mv_mult(pl_vec_out_out_int, pl_vec_in_in_int);

}
