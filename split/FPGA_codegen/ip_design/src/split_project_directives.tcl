# ################################################################## 
# Directives used by Vivado HLS project split_project

# DDR3 memory m_axi interface directives
set_directive_interface -mode m_axi -depth 35 "foo" memory_inout

# IP core handling directives
set_directive_interface -mode s_axilite -bundle BUS_A "foo"

# Input vectors offset s_axilite interface directives
set_directive_interface -mode s_axilite -register -bundle BUS_A "foo" byte_pl_vec_in_in_offset
set_directive_interface -mode s_axilite -register -bundle BUS_A "foo" byte_pl_mat_in_in_offset

# Output vectors offset s_axilite interface directives
set_directive_interface -mode s_axilite -register -bundle BUS_A "foo" byte_pl_vec_out_out_offset

set_directive_inline -off "foo_user"

# pipeline the for loop named "input_cast_loop_*" in foo.cpp
set_directive_pipeline "foo/input_cast_loop_pl_vec_in"
set_directive_pipeline "foo/input_cast_loop_pl_mat_in"
# pipeline the for loop named "output_cast_loop_*" in foo.cpp
set_directive_pipeline "foo/output_cast_loop_pl_vec_out"
