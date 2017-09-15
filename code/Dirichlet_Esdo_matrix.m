% Imposes (homogeneous) Dirichlet b.c. according to a list of nodes Node_flag
%  (on the matrix and rhs)
% Node_flag(k,1)=1 -> Dirichlet
% Node_flag(k,1)=2 -> Neumann
% !!!!! separate-displacement ordering !!!!!
%
function [S,K,A,B]=Dirichlet_Esdo_matrix(S,K,A,B,...
        Node_flagx,Node_flagy,nnode)
    
    nt  = size(S,1);
    flx = find(Node_flagx(:,1)==1);
    fly = find(Node_flagy(:,1)==1);
    dm  = mean(diag(S));
    
    S(flx,:) = 0;
    S(:,flx) = 0;
    S(fly+nnode,:) = 0;
    S(:,fly+nnode) = 0;
    
    K(flx,:) = 0;
    K(:,flx) = 0;
    K(fly+nnode,:) = 0;
    K(:,fly+nnode) = 0;
    
    A(flx,:) = 0;
    A(:,flx) = 0;
    A(fly+nnode,:) = 0;
    A(:,fly+nnode) = 0;
    
    B(flx,:) = 0;
    B(fly+nnode,:) = 0;
    
    dd = spalloc(nt,1,length(flx)+length(fly));
    dd(flx)=dm;
    dd(fly+nnode)=dm;
    
    S = S + spdiags(dd,0,nt,nt);
    K = K + spdiags(dd(1:2*nnode,1),0,2*nnode,2*nnode);
    return
