% Imposes (homogeneous) Dirichlet b.c. according to a list of nodes Node_flag
%  (on the matrix and rhs)
% Node_flag(k,1)=1 -> Dirichlet
% Node_flag(k,1)=2 -> Neumann
% !!!!! separate-displacement ordering !!!!!
%
function [rhs]=Dirichlet_Esdo_rhs(rhs,Node_flagx,Node_flagy,nnode,nnodeP)

rhs_d = rhs(1:2*nnode,1);
rhs_p = rhs(2*nnode+1:2*nnode+nnodeP);

flx = find(Node_flagx(:,1)==1);
fly = find(Node_flagy(:,1)==1);

rhs_d(flx,1)       = 0;
rhs_d(fly+nnode,1) = 0;
rhs = [rhs_d;rhs_p];

return
