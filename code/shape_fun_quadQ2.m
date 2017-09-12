% --------------------------------------------------------------------
% Evaluates the biquadratic b.f. at point 'k' from array 'points'
% 'points' contains the coordinates of the reference element
% (-1,-1),(1,1)
%
% Presummed ordering of the nodepoints in the reference f.e.:
% (-1,1)      (1,1)
%
%    2       6      3
%     ----------------      
%     |	       |        |	 
%     |        |        |
%   5|____9____|7      FUN(9)
%     |	       |       |	
%     |	       |       |  
%   1|____|____|4
%             8
%   
% (-1,-1)     (1,-1)
%
% f1d (g1d) are the 1D quadratic b.f. on the interval [-1, 0, 1]
%                                     point numbering [ 1  3  2]
% --------------------------------------------------------------------


function [BF]=shape_fun_quadQ2(ksi,eta)

%       ksi  = Gauss_point(1,k);
%       eta  = Gauss_point(2,k);
      
      th1 = inline('0.5*t.*(t-1)');
      th2 = inline('0.5*t.*(t+1)');
      th3 = inline('1-t.^2'); 
      
      BF(1,1) = th1(ksi).*th1(eta);
      BF(1,2) = th1(ksi).*th2(eta);
      BF(1,3) = th2(ksi).*th2(eta);
      BF(1,4) = th2(ksi).*th1(eta);
      BF(1,5) = th1(ksi).*th3(eta);
      BF(1,6) = th3(ksi).*th2(eta);
      BF(1,7) = th2(ksi).*th3(eta);
      BF(1,8) = th3(ksi).*th1(eta);
      BF(1,9) = th3(ksi).*th3(eta);

%       f1d(1) = 0.5*ksi*(ksi-1);
%       f1d(2) = 0.5*ksi*(ksi+1);
%       f1d(3) = 1-ksi*ksi;
%       g1d(1) = 0.5*eta*(eta-1);
%       g1d(2) = 0.5*eta*(eta+1);
%       g1d(3) = 1-eta*eta;
% 
%       vk=[1 1 2 2 1 3 2 3 3];
%       vl=[1 2 2 1 3 2 3 1 3];
%       for ii=1:9
%           k=vk(ii);
%           l=vl(ii);
%           BF(1,ii) = f1d(k)*g1d(l);
%       end
% %      BF(1,1:9) = f1d(vk).*g1d(vl);

return
