% Assembly of the element stiffness matrix for the 'rot' operator
% 
%
% S(18x18) = [ S11(9x9) S12(9x9) ]
%           [ S21(9x9) S22(9x9) ]
%
% G is assumed to have the following structure

%     [df_1/dx  df_2/dx df_3/dx ...df_9/dx]
%   G=[df_1/dy  df_2/dy df_3/dy ...df_9/dy]
%     [    *        *       *   ]

function [S]=assm_rot_quad(G)
    S11 = -G(2,:)'*G(2,:);
    S12 =  G(2,:)'*G(1,:); 
    S21 =  G(1,:)'*G(2,:);
    S22 = -G(1,:)'*G(1,:);
    S = zeros(18,18);
    S = [S11 S12;S21 S22];
return
% pointwise
%  for k=1:9,
%      for l=1:9,
% 	 S0(k,l)     = -G(2,k)*G(2,l);  %-dfk_y*dfl_y
% 	 S0(k,l+9)   = +G(2,k)*G(1,l);  %+dfk_y*dfl_x
% 	 S0(k+9,l)   = +G(1,k)*G(2,l);  %+dfl_x*dfk_y
% 	 S0(k+9,l+9) = -G(1,k)*G(1,l);  %-dfk_x*dfl_x
%      end
%  end
% 
%  
%  for k=1:9,
% 	 S1(k,1:9)     = -G(2,k)*G(2,1:9);  %-dfk_y*dfl_y
% 	 S1(k,10:18)   = +G(2,k)*G(1,1:9);  %+dfk_y*dfl_x
% 	 S1(k+9,1:9)   = +G(1,k)*G(2,1:9);  %+dfl_x*dfk_y
% 	 S1(k+9,10:18) = -G(1,k)*G(1,1:9);  %-dfk_x*dfl_x
%  end
