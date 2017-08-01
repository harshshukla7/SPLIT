classdef splitFuncTypes
 properties (Constant)
  vectorNorms = [splitFuncTypes.norm_one;...
   splitFuncTypes.norm_two;...
   splitFuncTypes.norm_inf];
  matrixNorms = [splitFuncTypes.norm_matrix_one;...
   splitFuncTypes.norm_matrix_two;...
   splitFuncTypes.norm_matrix_inf;...
   splitFuncTypes.norm_matrix_fro];
 end
 
 properties
  vexity
 end
 
 enumeration
  zero      (0)
  affine    (0)
  max       (1)
  min       (-1)
  abs       (1)
  quadratic (1)
  norm_one  (1)
  norm_two  (1)
  norm_inf  (1)
  norm_matrix_one  (1)
  norm_matrix_two  (1)
  norm_matrix_inf  (1)
  norm_matrix_fro  (1)
 end
 
 methods
  function type = splitFuncTypes(vexity_)
   type.vexity = vexity_;
  end
 end
 
 methods (Static)
  function str = toString(t)
   str = 'unkonwn';
   switch t
    case splitFuncTypes.zero
     str = 'zero';
    case splitFuncTypes.affine
     str = 'affine';
    case splitFuncTypes.max
     str = 'max';
    case splitFuncTypes.min
     str = 'min';
    case splitFuncTypes.abs
     str = 'absolute value';
    case splitFuncTypes.quadratic
     str = 'quadratic';
    case splitFuncTypes.norm_one
     str = 'one norm';
    case splitFuncTypes.norm_two
     str = 'two norm';
    case splitFuncTypes.norm_inf
     str = 'infinity norm';
    case splitFuncTypes.norm_matrix_one
     str = 'matrix one norm';
    case splitFuncTypes.norm_matrix_two
     str = 'matrix two norm';
    case splitFuncTypes.norm_matrix_inf
     str = 'matrix infinity norm';
    case splitFuncTypes.norm_matrix_fro
     str = 'frobenius norm';
   end
  end
 end
end
