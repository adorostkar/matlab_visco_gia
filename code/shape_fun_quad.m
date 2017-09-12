% --------------------------------------------------------------------
% Evaluates the bilinear b.f. at point 'k' from array 'points'
% 'points' contains the coordinates of the reference triangle
% (0,0),(1,0),(0,1)
%
% Presummed ordering of the nodepoints in the reference f.e.:
%  (-1,1) (1,1)
%     2    3
%      ----           
%     |    |         FUN(4)
%     |    |
%    1|____|4
%  (-1,-1)  (1,-1)
% --------------------------------------------------------------------


function [BF]=shape_fun_quad(ksi,eta)

%       ksi  = Gauss_point(1,k);
%       eta  = Gauss_point(2,k);
      ksim = (1-ksi);
      ksip = (1+ksi);
      etam = (1-eta);
      etap = (1+eta);

      BF(1,1) = 0.25*ksim*etam;
      BF(1,2) = 0.25*ksim*etap;
      BF(1,3) = 0.25*ksip*etap;
      BF(1,4) = 0.25*ksip*etam;

return
