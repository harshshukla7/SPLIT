

// Define FLOAT_FIX_VECTOR_NAME=1 to enable  fixed-point (up to 32 bits word length) arithmetic precision or 
// FLOAT_FIX_VECTOR_NAME=0 to enable floating-point single arithmetic precision.
#define FLOAT_FIX_STATE0_IN 0
#define FLOAT_FIX_PRIMAL0_IN 0
#define FLOAT_FIX_DUAL0_IN 0
#define FLOAT_FIX_TOL_ITERATES_IN 0
#define FLOAT_FIX_PRIMAL_OUT 0
#define FLOAT_FIX_DUAL_OUT 0
#define FLOAT_FIX_AUX_PRIMAL_OUT 0
#define FLOAT_FIX_ITERATES_OUT 0

//Input vectors INTEGERLENGTH:
#define STATE0_IN_INTEGERLENGTH 0
#define PRIMAL0_IN_INTEGERLENGTH 0
#define DUAL0_IN_INTEGERLENGTH 0
#define TOL_ITERATES_IN_INTEGERLENGTH 0
//Output vectors INTEGERLENGTH:
#define PRIMAL_OUT_INTEGERLENGTH 0
#define DUAL_OUT_INTEGERLENGTH 0
#define AUX_PRIMAL_OUT_INTEGERLENGTH 0
#define ITERATES_OUT_INTEGERLENGTH 0


//Input vectors FRACTIONLENGTH:
#define STATE0_IN_FRACTIONLENGTH 0
#define PRIMAL0_IN_FRACTIONLENGTH 0
#define DUAL0_IN_FRACTIONLENGTH 0
#define TOL_ITERATES_IN_FRACTIONLENGTH 0
//Output vectors FRACTIONLENGTH:
#define PRIMAL_OUT_FRACTIONLENGTH 0
#define DUAL_OUT_FRACTIONLENGTH 0
#define AUX_PRIMAL_OUT_FRACTIONLENGTH 0
#define ITERATES_OUT_FRACTIONLENGTH 0


//Input vectors size:
#define STATE0_IN_LENGTH 4
#define PRIMAL0_IN_LENGTH 46
#define DUAL0_IN_LENGTH 28
#define TOL_ITERATES_IN_LENGTH 4
//Output vectors size:
#define PRIMAL_OUT_LENGTH 46
#define DUAL_OUT_LENGTH 28
#define AUX_PRIMAL_OUT_LENGTH 28
#define ITERATES_OUT_LENGTH 4




typedef float data_t_memory;


#if FLOAT_FIX_STATE0_IN == 1
	typedef ap_fixed<STATE0_IN_INTEGERLENGTH+STATE0_IN_FRACTIONLENGTH,STATE0_IN_INTEGERLENGTH,AP_TRN,AP_WRAP> data_t_state0_in;
	typedef ap_fixed<32,32-STATE0_IN_FRACTIONLENGTH,AP_TRN,AP_WRAP> data_t_interface_state0_in;
#endif
#if FLOAT_FIX_PRIMAL0_IN == 1
	typedef ap_fixed<PRIMAL0_IN_INTEGERLENGTH+PRIMAL0_IN_FRACTIONLENGTH,PRIMAL0_IN_INTEGERLENGTH,AP_TRN,AP_WRAP> data_t_primal0_in;
	typedef ap_fixed<32,32-PRIMAL0_IN_FRACTIONLENGTH,AP_TRN,AP_WRAP> data_t_interface_primal0_in;
#endif
#if FLOAT_FIX_DUAL0_IN == 1
	typedef ap_fixed<DUAL0_IN_INTEGERLENGTH+DUAL0_IN_FRACTIONLENGTH,DUAL0_IN_INTEGERLENGTH,AP_TRN,AP_WRAP> data_t_dual0_in;
	typedef ap_fixed<32,32-DUAL0_IN_FRACTIONLENGTH,AP_TRN,AP_WRAP> data_t_interface_dual0_in;
#endif
#if FLOAT_FIX_TOL_ITERATES_IN == 1
	typedef ap_fixed<TOL_ITERATES_IN_INTEGERLENGTH+TOL_ITERATES_IN_FRACTIONLENGTH,TOL_ITERATES_IN_INTEGERLENGTH,AP_TRN,AP_WRAP> data_t_tol_iterates_in;
	typedef ap_fixed<32,32-TOL_ITERATES_IN_FRACTIONLENGTH,AP_TRN,AP_WRAP> data_t_interface_tol_iterates_in;
#endif
#if FLOAT_FIX_STATE0_IN == 0
	typedef float data_t_state0_in;
	typedef float data_t_interface_state0_in;
#endif
#if FLOAT_FIX_PRIMAL0_IN == 0
	typedef float data_t_primal0_in;
	typedef float data_t_interface_primal0_in;
#endif
#if FLOAT_FIX_DUAL0_IN == 0
	typedef float data_t_dual0_in;
	typedef float data_t_interface_dual0_in;
#endif
#if FLOAT_FIX_TOL_ITERATES_IN == 0
	typedef float data_t_tol_iterates_in;
	typedef float data_t_interface_tol_iterates_in;
#endif
#if FLOAT_FIX_PRIMAL_OUT == 1 
	typedef ap_fixed<PRIMAL_OUT_INTEGERLENGTH+PRIMAL_OUT_FRACTIONLENGTH,PRIMAL_OUT_INTEGERLENGTH,AP_TRN,AP_WRAP> data_t_primal_out;
	typedef ap_fixed<32,32-PRIMAL_OUT_FRACTIONLENGTH,AP_TRN,AP_WRAP> data_t_interface_primal_out;
#endif
#if FLOAT_FIX_DUAL_OUT == 1 
	typedef ap_fixed<DUAL_OUT_INTEGERLENGTH+DUAL_OUT_FRACTIONLENGTH,DUAL_OUT_INTEGERLENGTH,AP_TRN,AP_WRAP> data_t_dual_out;
	typedef ap_fixed<32,32-DUAL_OUT_FRACTIONLENGTH,AP_TRN,AP_WRAP> data_t_interface_dual_out;
#endif
#if FLOAT_FIX_AUX_PRIMAL_OUT == 1 
	typedef ap_fixed<AUX_PRIMAL_OUT_INTEGERLENGTH+AUX_PRIMAL_OUT_FRACTIONLENGTH,AUX_PRIMAL_OUT_INTEGERLENGTH,AP_TRN,AP_WRAP> data_t_aux_primal_out;
	typedef ap_fixed<32,32-AUX_PRIMAL_OUT_FRACTIONLENGTH,AP_TRN,AP_WRAP> data_t_interface_aux_primal_out;
#endif
#if FLOAT_FIX_ITERATES_OUT == 1 
	typedef ap_fixed<ITERATES_OUT_INTEGERLENGTH+ITERATES_OUT_FRACTIONLENGTH,ITERATES_OUT_INTEGERLENGTH,AP_TRN,AP_WRAP> data_t_iterates_out;
	typedef ap_fixed<32,32-ITERATES_OUT_FRACTIONLENGTH,AP_TRN,AP_WRAP> data_t_interface_iterates_out;
#endif
#if FLOAT_FIX_PRIMAL_OUT == 0 
	typedef float data_t_primal_out;
	typedef float data_t_interface_primal_out;
#endif
#if FLOAT_FIX_DUAL_OUT == 0 
	typedef float data_t_dual_out;
	typedef float data_t_interface_dual_out;
#endif
#if FLOAT_FIX_AUX_PRIMAL_OUT == 0 
	typedef float data_t_aux_primal_out;
	typedef float data_t_interface_aux_primal_out;
#endif
#if FLOAT_FIX_ITERATES_OUT == 0 
	typedef float data_t_iterates_out;
	typedef float data_t_interface_iterates_out;
#endif
