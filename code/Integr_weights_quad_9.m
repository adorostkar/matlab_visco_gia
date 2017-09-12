% --------------------------------------------------------------------
% Integr_weights_quad:
% Numerical integration using Gauss quadrature formulas
% for quadrilaterals - reference element data
% 9 integration points (higher accuracy)
% --------------------------------------------------------------------
% Presummed ordering of the nodepoints in the reference f.e.:
%  (-1,1) (1,1)
%     2    3
%      ----           
%     |    |       
%     |    |
%    1|____|4
%  (-1,-1)  (1,-1)
% --------------------------------------------------------------------

function [Gauss_point,Gauss_weight,np]=Integr_weights_quad_9

% Number of Gauss points
np = 9;
% Coordinates of Gauss points
vv=sqrt(3/5);
% Weights in Gauss points
w1 = 40/81;  %5/9*8/9
w2 = 25/81;  %5/9*5/9
w3 = 64/81;  %8/9*8/9

Gauss_point(1,1) = -vv; Gauss_point(2,1) = -vv;  Gauss_weight(1,1) = w2;
Gauss_point(1,2) =   0; Gauss_point(2,2) = -vv;   Gauss_weight(1,2) = w1;
Gauss_point(1,3) =  vv; Gauss_point(2,3) = -vv;  Gauss_weight(1,3) = w2;
Gauss_point(1,4) = -vv; Gauss_point(2,4) =   0;   Gauss_weight(1,4) = w1;
Gauss_point(1,5) =   0; Gauss_point(2,5) =   0;    Gauss_weight(1,5) = w3;
Gauss_point(1,6) =  vv; Gauss_point(2,6) =   0;   Gauss_weight(1,6) = w1;
Gauss_point(1,7) = -vv; Gauss_point(2,7) =  vv;  Gauss_weight(1,7) = w2;
Gauss_point(1,8) =   0; Gauss_point(2,8) =  vv;   Gauss_weight(1,8) = w1;
Gauss_point(1,9) =  vv; Gauss_point(2,9) =  vv;  Gauss_weight(1,9) = w2;


return
