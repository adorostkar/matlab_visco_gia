% Compute the element stiffness matrices corresponding to the
% elasticity problem in saddle-point form
% FUND - shape functions for the displacements
% DERD - derivatives of FUND
% FUNP - shape functions for the pressure
% DERP - derivatives of FUNP
% -----> in this particular case they are the same
% vec_coeff - ???? add description
% --------------------------------------------------------------------
function [E_elem,A_elem,M_elem9,L_elem9,R_elem9,G_elem9]=...
         Assm_ElAd_quadr_Q2(Gauss_point,Gauss_weight,...
                             FUNP_all,DERP_all,FUND_all,DERD_all,...
                             CoordP,CoordD,vec_coeff,wh)
                         
% global test_problem
global lan mju nju
      np    = 4;                    % number of points per f.e.
      npQ2  = 9;                    % number of points per f.e.
      nip   = size(Gauss_point,2);  % nip = number of integration points
      dim   = 2;                    % problem dimension
      ndof  = npQ2*dim;             % ndof = np*dim
      L_elem  = zeros(ndof,ndof);   % Laplace part
      A0_elem_b = zeros(ndof,ndof);  % advection of pre-stress
      A0_elem_c = zeros(ndof,ndof);  % buoyancy
      R_elem  = zeros(ndof,ndof);   % Rot part
      G_elem  = zeros(ndof,ndof);   % Grad-Div part
      M_elem  = zeros(npQ2,npQ2);   % Mass matrix for the elasticity (one diag block)
      E_elem  = zeros(ndof,ndof);   % Elasticity part
      E_elem1= zeros(ndof,ndof);   % Elasticity part

      A_elem  = zeros(ndof,ndof);   % Advection part 
      R_elem9 = zeros(ndof,ndof);   % Rot part only
      G_elem9 = zeros(ndof,ndof);   % Grad-Div part
      L_elem9 = zeros(ndof,ndof);   % displacement Laplace only
      M_elem9 = zeros(npQ2,npQ2);   % mass matrix displacements
      Z       = zeros(npQ2,npQ2);   % Zero block
% ------------------------------------------------------------------
% Gauss_point(2,nip)        Coord(nip,2)
% FUN(np)                   DER(2,np)          Jac(2,2)     IJac(2,2)
% Deriv(2,np)               B_cur(np,ndof)     D(np,np)
%
      for k=1:nip                         
         FUNP = FUNP_all{k};
         DERP = DERP_all{k};
         Jac   = DERP*CoordP;
	 IJac  = Jac\eye(2);                      % inv(JacP);
     Det   = det(Jac);
%         DerivP  = IJac*DERP;                      % (2x4)=(2x2)*(2x4)
%         M_elem4 = FUNP'*FUNP;                     % (4,4)=(4,1)*(1,4)
%         M0_elem = M0_elem + Det*Gauss_weight(k)*M_elem4;
%         L_elemc = DerivP'*DerivP;                 % (4,4)=(4,1)*(1,4)
%         L_elem4 = L_elem4 + Det*Gauss_weight(k)*L_elemc;

         FUND   = FUND_all{k};     % (9x1)
         FUNDT  = FUND';
         DERD   = DERD_all{k};     % (2x9)
         DerivD = IJac*DERD;                      % (2xnpQ2)=(2x2)*(2xnpQ2)
         L_elem0= DerivD'*DerivD; 
         L_elem = [L_elem0,Z;Z,L_elem0];
         
         G1_elem= DerivD(1,:)'*FUNP;              % (npQ2,np)=(npQ2,1)*(1,np)
         G2_elem= DerivD(2,:)'*FUNP;              % (npQ2,np)=(npQ2,1)*(1,np)
                  
         R_elem = assm_rot_quad(DerivD); 
	     G_elem = assm_graddiv_quad(DerivD);
         L_elem9= L_elem9 + Det*Gauss_weight(k)*L_elem; % Laplace terms only
         R_elem9= R_elem9 + Det*Gauss_weight(k)*R_elem; % rot terms only
         G_elem9= G_elem9 + Det*Gauss_weight(k)*G_elem; % grad-div terms only
         
%         E_elem = E_elem  + Det*Gauss_weight(k)*L_elem;% ... % (ndof x ndof)
         E_elem = E_elem...
                          +2*mju*Det*Gauss_weight(k)*L_elem... % (ndof x ndof)
                          +  mju*Det*Gauss_weight(k)*R_elem...  % rot terms
			                +  lan*Det*Gauss_weight(k)*G_elem;    % grad-div terms
%          E_elem1 = E_elem1...
%                           +mju*Det*Gauss_weight(k)*L_elem... % (ndof x ndof)
% 			                +  (mju+lan)*Det*Gauss_weight(k)*G_elem;    % grad-div terms
% 		if norm(E_elem1-E_elem)>1e-14, desktop, end	  
         A0_elem_b= advect_b(FUND,FUNDT,DerivD,Gauss_point(:,k),npQ2,vec_coeff,nju,wh); %pre-stress
         A0_elem_c= advect_c(FUND,FUNDT,DerivD,Gauss_point(:,k),npQ2,vec_coeff,nju,wh); %buoyancy
         A_elem = A_elem  + Det*Gauss_weight(k)*(-A0_elem_b + A0_elem_c); % advection terms                           
         M_elem  = FUNDT*FUND;                     % (9,9)=(9,1)*(1,9)
         M_elem9= M_elem9 + Det*Gauss_weight(k)*M_elem;
      end
      
return


% % The reference element
% CoordP=[-1 -1;-1 1;1 1;1 -1];
% CoordD=[-1 -1;-1 1;1 1;1 -1;-1 0;0 1; 1 0; 0 -1; 0 0];
% vec_coeff=[1,1,1,1];

% symbolic computations   
	b1 = vec_coeff(1);
	b2 = vec_coeff(2);
	c1 = vec_coeff(3);
	c2 = vec_coeff(4);
[MsymD,MsymP,LsymD,B1sym,B2sym,...
 R11sym,R12sym,...
 R21sym,R22sym,...
 G11sym,G12sym,...
 G21sym,G22sym,...
 Ad11sym,Ad12sym,...
 Ad21sym,Ad22sym,bsym]=THQ2Q1_symbolic(CoordP,CoordD,Gauss_point,b1,b2,c1,c2);   
    
%     
% 	[double(MsymP)-M0_elem]
    [double(MsymD)-M_elem9]
	[4*mju*double(LsymD)-L_elem9(1:9,1:9)]
% 	[double(B1sym)-B1_elem]
% 	[double(B2sym)-B2_elem]

	[double(Ad11sym) - A_elem(1:9,1:9)]
	[double(Ad12sym)  - A_elem(1:9,10:18)]
	[double(Ad21sym) - A_elem(10:18,1:9)]
	[double(Ad22sym) - A_elem(10:18,10:18)]          
	
% 	[double(R11sym)-R_elem9(1:9,1:9)]
% 	[double(R12sym)-R_elem9(1:9,10:18)]
% 	[double(R21sym)-R_elem9(10:18,1:9)]
% 	[double(R22sym)-R_elem9(10:18,10:18)]
%     
% 	[double(G11sym)-G_elem9(1:9,1:9)]
% 	[double(G12sym)-G_elem9(1:9,10:18)]
% 	[double(G21sym)-G_elem9(10:18,1:9)]
% 	[double(G22sym)-G_elem9(10:18,10:18)]

    [double(bsym)]
% 	
% %      APelem= advect_modified(A_elem);
% % subs(fiPdy,{x,y},{[-1/sqrt(3);-1/sqrt(3)],[1/sqrt(3);1/sqrt(3)]})	 

