function [Pos, Stress_cur] = SStime_midpoint(U,V,P,Node,Edge,Face,Face_eorder9,Face_Parent,nnode,hx,hy,time,Maxwell_time_inv)
% Compute approximate stress in the middle (center of mass, point 9) of each element.
% Use the derivatives of the basis functions.
% Assumed that the points in each face are ordered as in the reference
% element, thus we know explicitly the x- and y-orientation
% Strictly suited for rectandular meshes
% h is the fine mesh size
% Orientation of the face:
% (-1, 1)     (1, 1)
%
%     2       6       3
%     -----------------
%     |	        |        |
%     |         |        |
%    5|____9____|7
%     |	        |       |
%     |	        |       |
%    1|____|____|4
%              8
%
% (-1,-1)     (1,-1)
% -------------------------------

global E nju L_char M_char

Estar  = E/((1+nju)*(1-2*nju));
Coeff  = [1-nju nju 0; nju 1-nju 0; 0 0 (1-2*nju)/2];

mju = E/(2*(1+nju));
% U = L_char*U;
% V = L_char*V;
% Estar = M_char*Estar;

nface = size(Face_eorder9,2);
%hxi = 1/hx;
%hyi = 1/max(abs(hy));
 DERP_9 = shape_der_quad(0,0);    % DERP(2x4)
 DERD_9 = shape_der_quadQ2(0,0);  % DERV(2x9)
Pos = [];
for iface = 1:nface
    nl = Face_eorder9(:,iface);
    CoordP(1:4,1) = Node(1,nl(1:4))';  % Coord4(4,2)
    CoordP(1:4,2) = Node(2,nl(1:4))';
    CoordD(1:9,1) = Node(1,nl)';  % Coord9(9,2)
    CoordD(1:9,2) = Node(2,nl)';
    Jac   = DERP_9*CoordP;
    IJac  = Jac\eye(2);                      % inv(JacP);
    Det   = det(Jac);
    DerivD = IJac*DERD_9;
% x-derivative of the b.f.: DerivD(1)
% y-derivative of the b.f.: DerivD(2)
    dux = DerivD(1,:)*U(nl);
    duy = DerivD(2,:)*U(nl);
    dvx = DerivD(1,:)*V(nl);
    dvy = DerivD(2,:)*V(nl);

% Compute strain in vector form per node in the current face
% Stress_app(1:3,9)
    Strain(1,1) = dux;
    Strain(2,1) = dvy;
    Strain(3,1) = 0.5*(duy+dvx); %( with a factor 2 included)
    Stress = mean(P(nl(1:4)))*mju*[1; 1; 0] + 2*mju*Strain;
%     Stress =  Stress*(2-exp(-Maxwell_time_inv*time));
%     Stress(:) =  Estar*Coeff*Strain*(2-exp(-Maxwell_time_inv*time));
    Sx(iface) = Stress(1);
    Sy(iface) = Stress(2);
    Sxy(iface)= Stress(3);

    Pos = [Pos, Node(:, nl(9))];
end

Stress_cur = [Sx' Sy' Sxy'];

% for iface = 1:nface
%     nl = Face_eorder9(:,iface);
%     CoordD(1:9,1) = Node(1,nl)';  % Coord9(9,2)
%     CoordD(1:9,2) = Node(2,nl)';       
%     for k=1:9
%         A(k,:)=[1 CoordD(k,1) CoordD(k,2) CoordD(k,1)*CoordD(k,2) ...
% 	        CoordD(k,1)^2 CoordD(k,2)^2 CoordD(k,1)^2*CoordD(k,2) ...
% 	       CoordD(k,1)*CoordD(k,2)^2 CoordD(k,1)^2*CoordD(k,2)^2];
%     end 
%     A=A*1e4;
%     CO = A\eye(9);
%     CO = CO*1e4;
%     x9 = CoordD(9,1);
%     y9 = CoordD(9,2);
%     for k=1:9,
%         dfix(k) = CO(2,k) + CO(4,k)*y9 + 2*CO(5,k)*x9 + 2*CO(7,k)*x9*y9 + CO(8,k)*y9^2 + 2*CO(9,k)*x9*y9^2;
%         dfiy(k) = CO(3,k) + CO(4,k)*x9 + 2*CO(6,k)*y9 + CO(7,k)*x9^2 + 2*CO(8,k)*x9*y9 + 2*CO(9,k)*x9^2*y9;
%     end
%     dux = dfix*U(nl);
%     duy = dfiy*U(nl);
%     dvx = dfix*V(nl);
%     dvy = dfiy*V(nl);
%  % Compute strain in vector form per node in the current face
% % Stress_app(1:3,9)   
%     Strain1(1,1) = dux;
%     Strain1(2,1) = dvy;
%     Strain1(3,1) = duy+dvx; %( with a factor 2 included)
%     Stress1(:) =  Estar*Coeff*Strain1*(2-exp(-Maxwell_time_inv*time));
%     Sx1(iface) = Stress1(1);
%     Sy1(iface) = Stress1(2);
%     Sxy1(iface)= Stress1(3);
% end
% 
% Stress_cur1= [Sx1' Sy1' Sxy1'];

return
