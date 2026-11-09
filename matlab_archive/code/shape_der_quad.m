% --------------------------------------------------------------------
%
% Evaluates the first partial derivatives of the bilinear b.f.
% at point 'k' from array 'Gauss_point'
% 'Gauss_point' contains the coordinates of the
%  reference triangle (0,0),(1,0),(0,1)
%
% Presummed ordering of the nodepoints in the reference f.e.:
%  (-1,1) (1,1)
%     2    3
%      ----
%     |    |       DER = [ dfi_i/dksi; dfi_i/deta]
%     |    |       DER(2,4)
%    1|____|4
%  (-1,-1)  (1,-1)
% --------------------------------------------------------------------


function [DER]=shape_der_quad(ksi,eta) %
    
    %       ksi  = Gauss_point(1,k);
    %       eta  = Gauss_point(2,k);
    ksim = 1-ksi;
    ksip = 1+ksi;
    etam = 1-eta;
    etap = 1+eta;
    
    DER(1,1) = -0.25*etam;
    DER(2,1) = -0.25*ksim;
    
    DER(1,2) = -0.25*etap;
    DER(2,2) =  0.25*ksim;
    
    DER(1,3) =  0.25*etap;
    DER(2,3) =  0.25*ksip;
    
    DER(1,4) =  0.25*etam;
    DER(2,4) = -0.25*ksip;
    
    return
