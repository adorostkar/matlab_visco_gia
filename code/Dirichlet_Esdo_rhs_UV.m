% Imposes (homogeneous) Dirichlet b.c. according to a list of nodes Node_flag
%  (on the matrix and rhs)
% Node_flag(k,1)=1 -> Dirichlet
% Node_flag(k,1)=2 -> Neumann
% !!!!! separate-displacement ordering !!!!!
%
function [rhs_d]=Dirichlet_Esdo_rhs_UV(rhs,Node_flagx,Node_flagy,nnode)

rhs_d = rhs;

flx = find(Node_flagx(:,1)==1);
fly = find(Node_flagy(:,1)==1);

rhs_d(flx,1)       = 0;
rhs_d(fly+nnode,1) = 0;

return
