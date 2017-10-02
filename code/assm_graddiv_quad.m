
% Assembly of the element stiffness matrix for the 'grad div' operator
%  quadrilateral f.e., biquadratic b.f.
%
% S(18x18) = [ S11(9x9) S12(9x9) ]
%            [ S21(9x9) S22(9x9) ]
%
% G is assumed to have the following structure

%     [df_1/dx  df_2/dx df_3/dx ...df_9/dx]
%   G=[df_1/dy  df_2/dy df_3/dy ...df_9/dy]
%     [    *        *       *   ]

function [S]=assm_graddiv_quad(G)

           S11 = G(1,:)'*G(1,:);  %df_l/dx*df_k/dx

           S12 = G(1,:)'*G(2,:); %df_l/dx*df_k/dy

           S21 = G(2,:)'*G(1,:); %df_k/dx*df_l/dy

           S22 = G(2,:)'*G(2,:); %df_l/dy*df_k/dy

S = [S11 S12;
     S21 S22];

return
