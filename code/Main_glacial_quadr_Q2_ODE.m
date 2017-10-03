%     Standard formulation of the elasticity problem. Q2  discretization
% %%%%%%%%%%%%%%%%% 2D elasticity with advection terms %%%%%%%%%%%%%%%%%%%%%%%%%%
% Works on any given coarse mesh, described as {Node, Edge, Face}
% 1. The coarse mesh is refined 'levels' times
% 2. The mesh obtained in 1. is the pressure mesh. It is then refined
%    once more to get the mesh for the displacements. At this level
%    parent-child relation is kept to utilize the assembly of the
%    stiffness matrices
% 3. The viscoelastic solver is run. It includes the assembly of the
%    corresponding STIFFNESS and MASS matrices
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Block-factorized preconditioner, two-by-two block structure introduced
% Parent-child relation in the mesh refinement introduced
%
% ----- Description of the test problems
% test_problem = 0 ??
% test_problem = 8 ??
% test_problem = 9 Surface load on the whole top side
% test_problem = 10: GIA model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear Node* Edge* Face* Ice* uvp* UVP*
% % ------------------------------- global definitions
%global average_inner how_many outer_solv inner_solv inner_solver_type
%global restart restart_tol
%global avcgit_d avcgit_u
global h Maxwell_time_inv beta
global eps_inner
global theta
global rhs_const
global advec_const
global test_problem
global L_char M_char N_char U_char R_char G_char T_char T0

% Debug variable
% 1 - print debug texts
% 0 - no debug
global debug;
debug = 1;

% % Variables
% 0 - No additional output(figure and text)
% 1 - Only figures
% 2 - everything
verbose = 0;

np    = 4;  % number of points in the finite element
dim   = 2;  % number of degree of freedom (elasticity part)
dofP = 4;  % number of pressure degrees of freedom per element
dofD = 9;  % number of displacement degrees of freedom per element
ndof = dim*dofD;

% actionILU =[999,1e-1,1e-2,1e-3];
% actionILUt=['exact solve(A11)','cholinc(0.0100) ','cholinc(0.0010) ','cholinc(0.0001) '];
% actionFEM =[0,1,2,3,4,5,1e6,999];
% actionFEMt=['iter=0 ','iter=1 ','iter=2 ','iter=3 ','iter=4 ',...
%             'iter=5 ','iter=50','Z12    '];
% lll = 2;
% kkk = 1;

mkdir('./out')
diary('./out/diary')
% % Preparing parameters and mesh
tic
%% Test problem identifiers
% test_problem = 0; % pure elasticity in unit square, no surface  forces
% test_problem = 8; % visco-elasticity in unit square, no surface forces
test_problem = 9; % visco-elasticity in unit square, uniform  load 
% test_problem = 10; % visco-elasticity

switch test_problem
    case {0,8}
  wh = 'g0'; % manufactured solution
  beta = 0.1;   % exponent in the 'exact' solution
  no_domains = 2;
  Emagn = 1;
     [L,H,l_ice,h_ice,rho_ice,rho_earth,...
     Disco,Discoef,grav,load_surf,...
     T_LGM, T_EOG, T, delta_t_char] = Elasto_parameters(no_domains,wh,Emagn)
     H0=sign(H)*max(abs(H));
     Nx = 3;
     Ny = 3;
     [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_TH(L,H0,Nx,Ny);
   case 9
     wh = 'g9'; % uniform load on the top
      no_domains = 2;
      Emagn = 1; % can be 1, 10, 100 (jump in E between the two subdomains)
      beta = 0.1; % to be checked !!!
 % ------------------- unit square     
%      [L,H,l_ice,h_ice,rho_ice,rho_earth,...
%      Disco,Discoef,grav,load_surf,...
%      T_LGM, T_EOG, T, delta_t_char] = Elasto_parameters(no_domains,wh,Emagn);
%      H0=sign(H)*max(abs(H));
%       Nx = 3;
%       Ny = 3;

% ------------------- GIA geometry    test_case = 4;
          [L,H,l_ice,h_ice,rho_ice,rho_earth,...
          Disco,Discoef,grav,load_surf,...
          T_LGM, T_EOG, T, delta_t_char] = Visco_parameters_const_load(no_domains,wh,Emagn);
           H0=-max(abs(H));
          Nx = 11; %L/l_ice+1;      % ensure mesh aligned with the ice after one refinement
          Ny = 5; %abs(H0)/l_ice+1;
       [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_TH(L,H0,Nx,Ny);
    case 10
     wh = 'gs';  Maxwell_time_inv =0;
     no_domains = 2;
     Emagn = 1; % can be 1, 10, 100 (jump in E between the two subdomains)
     test_case = 4;    
    switch test_case
    case 1  % Bjorn
    [L,H,l_ice,h_ice,rho_ice,rho_earth,...
       Disco,Discoef,grav,load_surf,...
       T_LGM, T_EOG, T, delta_t_char] = Visco_parameters_Bjorn(no_domains,wh,Emagn);
       H0=-max(abs(H));
       Nx = 11;%L/(2*l_ice)+1;      % ensure mesh aligned with the ice after one refinement
       Ny = 5;%abs(H0)/(2*l_ice)+1;
       [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_TH_Bjorn(L,H0,Nx,Ny);
    case 2  % Wu 1992
       [L,H,l_ice,h_ice,rho_ice,rho_earth,...
          Disco,Discoef,grav,load_surf,...
          T_LGM, T_EOG, T, delta_t_char] = Visco_parameters_Wu92(no_domains,wh,Emagn);
       H0=-max(abs(H));
       Nx = 11;%L/(2*l_ice)+1;      % ensure mesh aligned with the ice after one refinement
       Ny = 5;%abs(H0)/(2*l_ice)+1;

       [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_TH_wu92(L,H0,Nx,Ny);
%  [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_Wu_coarse(L,H,l_ice);
    case 3  % Wu 1998
       [L,H,l_ice,h_ice,rho_ice,rho_earth,...
          Disco,Discoef,grav,load_surf,...
          T_LGM, T_EOG, T, delta_t_char] = Visco_parameters_Wu98(no_domains,wh,Emagn);
       H0=-max(abs(H));
       Nx = 51;  % ensure mesh aligned with the ice after one refinement
       Ny = 9;
       [xc,yc,hx,hx2,hy,Nx,Ny] = Glace_coord_vectors_TH_Wu98(L,H0,Nx,Ny);
%  [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_Wu_coarse(L,H,l_ice);
    otherwise  % standard testing
       [L,H,l_ice,h_ice,rho_ice,rho_earth,...
          Disco,Discoef,grav,load_surf,...
          T_LGM, T_EOG, T, delta_t_char] = Visco_parameters_const_load(no_domains,wh,Emagn);
       H0=-max(abs(H));
       Nx = L/l_ice+1;      % ensure mesh aligned with the ice after one refinement
       Ny = abs(H0)/l_ice+1;
       [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_TH(L,H0,Nx,Ny);
    end
end

[Node,Edge,Face,...
 Node_flagx,Node_flagy,...
 Edge_flagx,Edge_flagy,...
 Face_flag,Face_thick] = Rectan_glace_vect(L,H,xc,yc,Nx,Ny,...
                                           no_domains,Disco,wh);


% Visualise the mesh
if(verbose ~= 0)
    figure(1),clf,Bvisual_mesh(Node,Edge,Face,1,1,1,3,16)
end
disp(['Time to create initial mesh: ' num2str(toc)])

% -------------------- Input parameters ---------
nrefin = input('How many times to refine: ');
   hx = hx/2^nrefin;
   hy = hy/2^nrefin;
   h  = min(abs(hx),min(abs(hy)));
   sigma = h^2;
%    disp('sigma  = h^2')

solver_type = 1;
levels    = nrefin + 1;
lvl_total = levels + 1;
lvl_coars = max(1,levels-1);
nface_lvl = zeros(lvl_total,1);
nnode_lvl = zeros(lvl_total,1);
nface_lvl(1) = size(Face,2);
nnode_lvl(1) = size(Node,2);

tic
for lvl=1:nrefin,
    [Node,Edge,Face,...
     Node_flagx,Node_flagy,...
     Edge_flagx,Edge_flagy,...
     Face_flag,Face_thick] = my_Refine_quadr(Node,Edge,Face,...
                                   Node_flagx,Node_flagy,...
                                   Edge_flagx,Edge_flagy,...
                                   Face_flag,Face_thick);
     nface_lvl(lvl+1) = size(Face,2);
     nnode_lvl(lvl+1) = size(Node,2);
%     figure(1),Bvisual_mesh(Node,Edge,Face,1,1,1,0,16)
end

nnodeP = size(Node,2);
nedgeP = size(Edge,2);
nfaceP = size(Face,2);    % number of subdomains (for the pressure)
nallP  = nnodeP;          % number of pressure variables
% disp(['Total number of quadrilaterals for the pressure: ' int2str(nfaceP)])

% figure(1),clf,Bvisual_mesh(Node,Edge,Face,1,1,1,1,16)
disp(['Total number of coarse quadrilaterals: ' int2str(nfaceP)])

%% ----- detect boundary edges under the load (ice) based on the pressure mesh
xmax = max(Node(1,:));
xmin = min(Node(1,:));
ymax = max(Node(2,:));
ymin = min(Node(2,:));

% ---------------- Create mesh-related matrices (again) -----
nnode = nnodeP+nedgeP+nfaceP;
nedge = 2*nedgeP+4*nfaceP;
nface = 4*nfaceP;      % number of subdomains (for the displacements)
Edge_Node = spalloc(nedgeP,nnodeP,2*nedgeP);
Face_Edge = spalloc(nfaceP,nedgeP,4*nedgeP);
for iedge=1:nedgeP,
%     vf = Edge(:,iedge)';
%     Edge_Node(iedge,vf)=1;
    Edge_Node(iedge,Edge(1,iedge))=1;
    Edge_Node(iedge,Edge(2,iedge))=1;
end
for iface=1:nfaceP,
%     vf=Face(:,iface)';
%     Face_Edge(iface,vf)=1;
    Face_Edge(iface,Face(1,iface))=1;
    Face_Edge(iface,Face(2,iface))=1;
    Face_Edge(iface,Face(3,iface))=1;
    Face_Edge(iface,Face(4,iface))=1;
end

    Load_Nodes=[];  Load_Edges=[];  Load_Edges_list=[];  Load_Faces=[];
    Top_Nodes=[];  Top_Edges=[];  Top_Edges_list=[];  Top_Faces=[];
    Right_Nodes=[];Right_Edges=[];Right_Edges_list=[];Right_Faces=[];        


switch test_problem
    case{0,8}  % no surface load imposed on the boundary of the domain
        % nothing to do
    case{9,10} % load on the top boundary (uniform or partial)
        Surface_Nodes = find(Node(2,:)==ymax); % all surface nodes
% ~~~~ Find the part where the load is applied (start)
        [~,noj]   = find(Node(1,Surface_Nodes(:))<=l_ice);
        Load_Nodes = Surface_Nodes(noj)
        [noi,noj]=find(Edge_Node(:,Load_Nodes));
        noi=unique(noi);
        clear Load_edges, lb=0;
        for i=1:length(noi),
            if (Node(2,Edge(:,noi(i)))==[ymax ymax])&(prod(Node(1,Edge(:,noi(i))))<=l_ice^2),
                lb = lb + 1; Load_Edges_list(lb)=noi(i);
            end
        end
        Load_Edges=Edge(:,Load_Edges_list);
%----   Find faces under the ice
        Edge_Face = Face_Edge';
        [wv,Load_Faces] = find(Edge_Face(Load_Edges_list,:));
        Load_Faces(wv) = Load_Faces;
   % ~~~~ Find the part where the load is applied (end)
   % ~~~~ Find the top boundary  no load (start)
        [~,noj]   = find(Node(1,Surface_Nodes(:))>l_ice);
        if ~isempty(noj)
            Top_Nodes = Surface_Nodes(noj);
            [noi,noj]=find(Edge_Node(:,Top_Nodes));
            noi=unique(noi);
           clear Top_edges, lb=0;
          for i=1:length(noi),
              if (Node(2,Edge(:,noi(i)))==[ymax ymax])&(prod(Node(1,Edge(:,noi(i))))>l_ice^2),
                  lb = lb + 1; Top_Edges_list(lb)=noi(i);
              end
          end
        Top_Edges=Edge(:,Top_Edges_list);
%----   Find faces under the ice-free surface
         [wv,Top_Faces] = find(Edge_Face(Top_Edges_list,:));
         Top_Faces(wv) = Top_Faces;
        end
    end % switch test_problem

% Refine once to obtain the mesh for the displacements and
% Face_eorder4 and Face_eorder9 arrays

   [Node,Edge,Face,...
    Node_flagx,Node_flagy,...
    Edge_flagx,Edge_flagy,...
    Face_flag,Face_thick,...
    Face_Node,Face_Parent,...
    Face_eorder9,Face_eorder4,...
    nface_lvl,nnode_lvl] = my_Refine_quadr_hier(Node,Edge,Face,...
                           Node_flagx,Node_flagy,Edge_flagx,Edge_flagy,...
                           Face_flag,Face_thick,nface_lvl,nnode_lvl,levels);
%  figure(1),clf,Bvisual_mesh(Node,Edge,Face,1,1,1,3,16)
   hx = hx/2;
   hy = hy/2;

nnode = size(Node,2);
nedge = size(Edge,2);
nface = size(Face,2);      % number of subdomains (for the displacements)

% ---------------- Update mesh-related matrices (again) -----
Edge_Node = spalloc(nedge,nnode,2*nedge);
Face_Edge = spalloc(nface,nedge,4*nedge);
for iedge=1:nedge,
%     vf = Edge(:,iedge)';
%     Edge_Node(iedge,vf)=1;
    Edge_Node(iedge,Edge(1,iedge))=1;
    Edge_Node(iedge,Edge(2,iedge))=1;
end
for iface=1:nface,
%     vf=Face(:,iface)';
%     Face_Edge(iface,vf)=1;
    Face_Edge(iedge,Face(1,iface))=1;
    Face_Edge(iedge,Face(2,iface))=1;
    Face_Edge(iedge,Face(3,iface))=1;
    Face_Edge(iedge,Face(4,iface))=1;
end
% ------------ update surface nodes
switch test_problem
    case 0
        surface_Nodes=[];
    case 8 % the surface load is on the bottom boundary
    Surface_Nodes = find(Node(2,:)==ymin); % all surface nodes
    case {9,10}
    Surface_Nodes = find(Node(2,:)==ymax); % all surface nodes
end
% ------------ detect boundary edges
%bndry_edge = zeros(nedge,1);
bndry_edge = sum(Face_Edge,1);    % The boundary edges are with sum '1'
[noi,Bndry_Edges]=find(bndry_edge==1); % Bndry_Edges is a list of boundary edges ONLY!

%----->  Clear unnecessary arrays
clear Face_Edge Face_Node Edge_Node

disp(['End refinement. Elapsed: ' num2str(toc)])
if size(Node,2)<5, figure(1),clf,Bvisual_mesh(Node,Edge,Face,0,0,0,3,12), end

disp(' Ready with the mesh.')

disp(' Start the visco-elastic solver.')
Tmax = T;%input(' Run till time: ');

eps_inner = 5e-1;
Visco_elastic_solver_Q2_ODE

return
