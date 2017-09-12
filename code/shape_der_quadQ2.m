% --------------------------------------------------------------------
%
% Evaluates the first partial derivatives of the biquadratic b.f.
% at point 'k' from array 'Gauss_point'
% 'Gauss_point' contains the coordinates of the
%  reference triangle (0,0),(1,0),(0,1)
%
% (-1,-1),(1,1)
%
% Presummed ordering of the nodepoints in the reference f.e.:
% (-1,1)      (1,1)
%
%    2    6     3
%     ---------
%    |	  |    |
%    |	  |    |
%   5|____9____|7      FUN(9)
%    |	  |    |
%    |	  |    |
%   1|____|____|4
%         8
%
% (-1,-1)     (1,-1)
% --------------------------------------------------------------------



function [DER]=shape_der_quadQ2(ksi,eta) %
    
    %       ksi  = Gauss_point(1,k);
    %       eta  = Gauss_point(2,k);
    ksim = (ksi-0.5);
    ksip = (ksi+0.5);
    etam = (eta-0.5);
    etap = (eta+0.5);
    
    DER(1,1) = 0.5*eta*(etam - 0.5)*ksim;
    DER(2,1) = 0.5*ksi*(ksim - 0.5)*etam;
    
    DER(1,2) = 0.5*eta*(etap + 0.5)*ksim;
    DER(2,2) = 0.5*ksi*(ksim - 0.5)*etap;
    
    DER(1,3) = 0.5*eta*(etap + 0.5)*ksip;
    DER(2,3) = 0.5*ksi*(ksip + 0.5)*etap;
    
    DER(1,4) = 0.5*eta*(etam - 0.5)*ksip;
    DER(2,4) = 0.5*ksi*(ksip + 0.5)*etam;
    
    DER(1,5) = -ksim*(eta^2 - 1);
    DER(2,5) = -ksi*eta*(ksim - 0.5);
    
    DER(1,6) = -ksi*eta*(etap + 0.5);
    DER(2,6) = -etap*(ksi^2 - 1);
    
    DER(1,7) = -ksip*(eta^2 - 1);
    DER(2,7) = -ksi*eta*(ksip + 0.5);
    
    DER(1,8) = -ksi*eta*(etam - 0.5);
    DER(2,8) = -etam*(ksi^2 - 1);
    
    DER(1,9) = 2*ksi*(eta^2 - 1);
    DER(2,9) = 2*eta*(ksi^2 - 1);
    
    return
    
    %  x=ksi; y=eta;
    %  df(1,1) = 0.5*y*(y - 1)*(x - 0.5);
    %  df(2,1) = 0.5*y*(y + 1)*(x - 0.5);
    %  df(3,1) = 0.5*y*(y + 1)*(x + 0.5);
    %  df(4,1) = 0.5*y*(y - 1)*(x + 0.5);
    % df(5,1) = (x - 0.5)*(1 - y^2 );
    % df(6,1) = -x*y*(y + 1);
    % df(7,1) = (x + 0.5)*(1 - y^2 );
    % df(8,1) = x*y*(y - 1);
    % df(9,1) = 2*x*(y^2 - 1);
    %
    % df(1,2) = -0.5*x*(1 - x)*(y - 0.5);
    % df(2,2) = -0.5*x*(1 - x)*(y + 0.5);
    % df(3,2) = 0.5*x*(x + 1)*(y + 0.5);
    % df(4,2) = 0.5*x*(x + 1)*(y - 0.5);
    % df(5,2) = y*x*(1 - x);
    % df(6,2) = -(y + 0.5)*(x^2 - 1);
    % df(7,2) = -y*x*(x + 1);
    % df(8,2) = -(y - 0.5)*(x^2 - 1);
    % df(9,2) = 2*y*(x^2 - 1);
    %
    % return
