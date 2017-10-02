% Given: Node,Edge,Face
% asembling of the global stiffness matrix (Neumann everywhere)
% K(nnode,nnode) and the right-hand-side vector (only body forces)
% Standard formulation

function [K,A,L,M,rhs_d] = Assembly_ElAd_quadr_Q2(Node,...
                           Face_eorder9,Face_eorder4,...
                           Face_flag,Face_thick,Disco,Discoef,...
                           Gauss_point,Gauss_weight,...
                           FUNP_all,DERP_all,FUND_all,DERD_all,...
                           nnode_lvl,nface_lvl,lvl_total,vec_coeff,wh)
global verbose
global advec_const
global lan mju nju

nnodeD = nnode_lvl(lvl_total);   %number of nodes for the displacements
nnodeP = nnode_lvl(lvl_total-1); %number of nodes for the pressure
nface  = nface_lvl(lvl_total-1); %number of faces for the pressure
dim   = 2;	   %
dofd =  9; 
dofD = dofd*dim; % ndof = nip*dim
dofP =  4;

nrmmax = 0;
K  = spalloc(nnodeD,nnodeD,9*nnodeD); % velocity stiffness matrix
A  = spalloc(nnodeD,nnodeD,9*nnodeD); % convection matrix
rhs_d= zeros(2*nnodeD,1);
nallV = 2*nnodeD;
nallP = nnodeP;


if(verbose ~= 0)
    disp('Begin allocating memory...')
end

lengthK=nface*dofD*dofD;
KI  = zeros(lengthK,1);
KJ  = zeros(lengthK,1);
KV  = zeros(lengthK,1);

LI  = zeros(lengthK,1);  % Laplace only
LJ  = zeros(lengthK,1);
LV  = zeros(lengthK,1);

RI  = zeros(lengthK,1);  % Rot only
RJ  = zeros(lengthK,1);
RV  = zeros(lengthK,1);

GI  = zeros(lengthK,1);  % Grad-Div only
GJ  = zeros(lengthK,1);
GV  = zeros(lengthK,1);

MI  = zeros(lengthK,1);  % Laplace only
MJ  = zeros(lengthK,1);
MV  = zeros(lengthK,1);
nextK = 0;

lengthA=nface*dofD*dofD;  % Pre-stress advection only
AI  = zeros(lengthA,1);
AJ  = zeros(lengthA,1);
AV  = zeros(lengthA,1);
nextA = 0;



if(verbose ~= 0)
    disp('...end allocating memory.')
end
			  
for iface_p=1:nface,     % ---> walk on the parent (pressure) elements
%    disp(['Pressure face ' int2str(iface_p)])
     local_node_listP = Face_eorder4(:,iface_p);
     local_node_listD0= Face_eorder9(:,iface_p);
     local_node_listD=order_face_nodes_quadQ2(Node,local_node_listP,local_node_listD0);
     if local_node_listD0~=local_node_listD, 
         pause
     end

   CoordP(1:4,1) = Node(1,local_node_listP)';  % Coord4(4,2)
   CoordP(1:4,2) = Node(2,local_node_listP)';
   CoordD(1:9,1) = Node(1,local_node_listD)';  % Coord9(9,2)
   CoordD(1:9,2) = Node(2,local_node_listD)';

% ------------------------ determine nju, E
nju = Discoef(1,Face_flag(iface_p,1));
E   = Discoef(2,Face_flag(iface_p,1))*Face_thick(iface_p,1);

% ------------------------ compute the coefficients lan, mju
mju = E/(2*(1+nju));        lan = 2*mju*nju/(1-2*nju);
%rho = (1-2*nju)/(2*nju);    % mju/lan;

% - - - - - generation of the element matrices:

    [El_elem,Ad_elem,M_elem,L_elem9,R_elem9,G_elem9]=...
     Assm_ElAd_quadr_Q2(Gauss_point,Gauss_weight,...
                             FUNP_all,DERP_all,FUND_all,DERD_all,...
                             CoordP,CoordD,vec_coeff,wh);

% -------------- the following ordering is NOT for separate displacements
% local_node_liste = [2*local_node_list(1)-1,2*local_node_list(1),...
%                     2*local_node_list(2)-1,2*local_node_list(2),...
%                     2*local_node_list(3)-1,2*local_node_list(3),...
%                     2*local_node_list(4)-1,2*local_node_list(4)];

% -------------- the following ordering is for separate displacements     
% local_node_liste = [local_node_list,local_node_list+nnode];	

% local_node_liste9= [macro_node_list(1),macro_node_list(2),...
%                     macro_node_list(3),macro_node_list(4),...
%                     macro_node_list(1)+9,macro_node_list(2)+9,...
%                     macro_node_list(3)+9,macro_node_list(4)+9];	


macro_node_listD=[local_node_listD;local_node_listD+nnodeD]; % sdo
ZZ = zeros(size(M_elem));
M_elem9 = [M_elem ZZ;ZZ M_elem];
if ~strcmp(wh,'g0')
% % No body forces in this case  
   rhs_elemD = zeros(dofD,1);
else
   bforce = bodyf_time(Node(1,local_node_listD)',Node(2,local_node_listD)',0,...
                  mju,lan,vec_coeff);          
% % - - - - - bforce is a (ndof x 1) vector (for quadrilaterals)
% %           IN separate displacements ordering !!!
     rhs_elemD = sum(M_elem9)'.*bforce;  % M_elem9*bforce;
end 

% Assembly of the global and macro-element stiffness and mass matrices
     for ii=1:dofD,
       is = macro_node_listD(ii);
       for jj=1:dofD,
           js = macro_node_listD(jj);
%          K(is,js) = K(is,js) + mju*El_elem(ii,jj);
           nextK = nextK + 1;
           KI(nextK) = is;
           KJ(nextK) = js;
           KV(nextK) = El_elem(ii,jj);

           LI(nextK) = is;
           LJ(nextK) = js;
           LV(nextK) = L_elem9(ii,jj);     % vector Laplace
           RI(nextK) = is;
           RJ(nextK) = js;
           RV(nextK) = R_elem9(ii,jj);     % rot
           GI(nextK) = is;
           GJ(nextK) = js;
           GV(nextK) = G_elem9(ii,jj);     % grad-div
           MI(nextK) = is;
           MJ(nextK) = js;
           MV(nextK) = M_elem9(ii,jj);            % mass matrix

           nextA = nextA + 1;
           AI(nextA) = is;
           AJ(nextA) = js;
           AV(nextA) = Ad_elem(ii,jj); %all advection coeffs must already 
                                       %be included in the element matrix
       end
       rhs_d(is,1) = rhs_d(is,1) + rhs_elemD(ii,1);
    end    
end             % end pressure iface_p-loop

KI = KI(1:nextK,1); KJ = KJ(1:nextK,1); KV = KV(1:nextK,1);
LI = LI(1:nextK,1); LJ = KJ(1:nextK,1); LV = LV(1:nextK,1);
MI = MI(1:nextK,1); MJ = MJ(1:nextK,1); MV = MV(1:nextK,1);
AI = AI(1:nextA,1); AJ = AJ(1:nextA,1); AV = AV(1:nextA,1); % 

K  = sparse(KI, KJ, KV);
L  = sparse(LI, LJ, LV);
M  = sparse(MI, MJ, MV);
A  = advec_const*sparse(AI, AJ, AV);

return



