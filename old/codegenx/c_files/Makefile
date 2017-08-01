CFLAGS = -O3 -Os

ama: solver_wrapper.c solver_ama.c basic_operations.c projections.c solve_KKT_ama.c compute_vector_to_project.c compute_over_relaxation.c compute_projection.c compute_primal_residual.c compute_dual_residual.c 
	gcc $(CFLAGS) -o ama solver_wrapper.c solver_ama.c basic_operations.c projections.c solve_KKT_ama.c compute_vector_to_project.c compute_over_relaxation.c compute_projection.c compute_primal_residual.c compute_dual_residual.c

admm: solver_wrapper.c solver_admm.c basic_operations.c projections.c solve_KKT_admm.c compute_vector_to_project.c compute_over_relaxation.c compute_projection.c compute_primal_residual.c compute_dual_residual.c
	gcc $(CFLAGS) -o admm solver_wrapper.c solver_admm.c basic_operations.c projections.c solve_KKT_admm.c compute_vector_to_project.c compute_over_relaxation.c compute_projection.c compute_primal_residual.c compute_dual_residual.c 

fama: solver_wrapper.c solver_fama.c basic_operations.c projections.c solve_KKT_fama.c compute_vector_to_project.c compute_over_relaxation.c compute_projection.c compute_primal_residual.c compute_dual_residual.c 
	gcc $(CFLAGS) -o fama solver_wrapper.c solver_fama.c basic_operations.c projections.c solve_KKT_fama.c compute_vector_to_project.c compute_over_relaxation.c compute_projection.c compute_primal_residual.c compute_dual_residual.c 

fadmm: solver_wrapper.c solver_fadmm.c basic_operations.c projections.c solve_KKT_fadmm.c compute_vector_to_project.c compute_over_relaxation.c compute_projection.c compute_primal_residual.c compute_dual_residual.c 
	gcc $(CFLAGS) -o fadmm solver_wrapper.c solver_fadmm.c basic_operations.c projections.c solve_KKT_fadmm.c compute_vector_to_project.c compute_over_relaxation.c compute_projection.c compute_primal_residual.c compute_dual_residual.c 

CPI: solver_wrapper.c solver_cp1.c basic_operations.c projections.c compute_array_p.c compute_array_nu.c compute_projection_conjugates.c compute_array_x.c compute_array_xbar.c compute_primal_residual_chambolle.c compute_dual_residual_chambolle.c 
	gcc $(CFLAGS) -o CPI solver_wrapper.c solver_cp1.c basic_operations.c projections.c compute_array_p.c compute_array_nu.c compute_projection_conjugates.c compute_array_x.c compute_array_xbar.c compute_primal_residual_chambolle.c compute_dual_residual_chambolle.c

CPII: solver_wrapper.c solver_cp2.c basic_operations.c projections.c compute_array_p.c compute_array_nu.c compute_projection_conjugates.c compute_array_x_varying_stepsize.c compute_array_xbar.c compute_primal_residual_chambolle.c compute_dual_residual_chambolle.c 
	gcc $(CFLAGS) -o CPII solver_wrapper.c solver_cp2.c basic_operations.c projections.c compute_array_p.c compute_array_nu.c compute_projection_conjugates.c compute_array_x_varying_stepsize.c compute_array_xbar.c compute_primal_residual_chambolle.c compute_dual_residual_chambolle.c
    
clean:
	rm -rf test_admm test_ama test_fama test_fadmm ama admm fama fadmm CPI CPII

