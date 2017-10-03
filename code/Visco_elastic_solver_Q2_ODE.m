%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The viscoelastic solver - ODE
%
%
%[Gauss_pointP,Gauss_weightP] = Intnallegr_weights_quad;
[Gauss_point,Gauss_weight,nipD] = Integr_weights_quad_9;

% get the order of the nodes
[X,Y, id] = reorder_nodes(Node);

FUND_all=cell(9,1);
DERD_all=cell(9,1);
for k=1:nipD
% 9-point integration f-la
     FUNP_all{k} = shape_fun_quad(Gauss_point(1,k),Gauss_point(2,k));    % FUNP(1,4)
     DERP_all{k} = shape_der_quad(Gauss_point(1,k),Gauss_point(2,k));    % DERP(2x4)
     FUND_all{k} = shape_fun_quadQ2(Gauss_point(1,k),Gauss_point(2,k));  % FUNV(1,9)
     DERD_all{k} = shape_der_quadQ2(Gauss_point(1,k),Gauss_point(2,k));  % DERV(2x9)
end

T_BEG = 0; % the time point from which the ise load is imposed

nallD = 2*nnode;
nall  = nallD;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocations
Face_estiffS = []; %zeros(ndof+np,ndof+np,nface);

S_char  = 1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E1. Assembly the FE matrix,  K,A,B,C,
%     S = [K+A B^T]
%         [ B  -C]
%      elasticity material coefficients at time tau=t
%      initial load f (body forces are zero)
% K - elasticity part of S11
% A - advection part of S11  (S11=K+A)
%     [K+A     [B1 B2]^T ]
% S = [[B1 B2]   -C      ]
% K and A - in separate displacement ordering
% C - pressure mass matix

Maxwell_time_inv = Discoef(3,1);%

if test_problem < 10,
%    disp('NO PRE-STRESS ADVECTION.')  % no advection
   vec_coeff = [0,0,0,0]
   disp('EXACT SOLUTION WITH ADVECTION.')  %
%    vec_coeff = [1,2,3,4]
else
   disp('WITH PRE-STRESS ADVECTION.')  % with advection
 vec_coeff = [0,0,0,0]
%   vec_coeff = [0,-rho_earth*grav,0,-rho_earth*grav]
%  vec_coeff = [0,-rho_earth*grav,0,0]
end
% vec_coeff = vec_coeff0.*advec_const  % scaling of the advection terms
k = 1;
secs_per_year = 365*24*3600; % Duration of a year in seconds 3.1536e+7
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% VERY UGLY FIX
global E nju mju lan 
    nju = Discoef(1,Face_flag(1,1));
    E   = Discoef(2,Face_flag(1,1))*Face_thick(1,1);
    mju = E/(2*(1+nju));        lan = 2*mju*nju/(1-2*nju);
 [K0,A0,L0,M,rhs_b] = Assembly_ElAd_quadr_Q2(Node,...
                           Face_eorder9,Face_eorder4,...
                           Face_flag,Face_thick,Disco,Discoef,...
                           Gauss_point,Gauss_weight,...
                           FUNP_all,DERP_all,FUND_all,DERD_all,...
                           nnode_lvl,nface_lvl,lvl_total,vec_coeff,wh); %pure elastic
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_cur0 = K0+A0;  % y-cordinate axis pointing upwards
if (verbose ~= 0)
   disp('PRE-STRESS ADVECTION MATRIX WITH ''+''.')  %
end

% rhs_exact0 = bodyf_time(Node(1,:)',Node(2,:)',0, mju,lan,vec_coeff);
% rhs_exact = M*rhs_exact0;
% rhs_exactl= sum(M)'.*rhs_exact0;

time_cur = 0;                      %  rhs & b.c. for the elastic response

%   Add surface forces
[rhs_s] = Assembly_rhs_Q2(Node,Edge,...
                               wh,l_ice,h_ice,rho_ice,grav,...
                               Load_Edges,Load_Edges_list,...
                               nnodeP,time_cur,T_BEG,T_LGM,T_EOG);
                           
%  Compute the initial rhs0 (body forces + surface tractions)
rhs_cur0 = rhs_b + rhs_s;

% Dirichlet_visco_UV           %  impose b.c. both to T_cur and rhs_cur
%% ------------------ impose Dirichlet b.c.; standard formulation ------------
if wh=='g0'    % ... exact solution known, homogeneous material assumed
    [Uex,Vex] = UVtime_sol(Node(1,:)',Node(2,:)',time_cur); 
    uv_sol = [Uex;Vex];
    [T_cur,rhs_cur] = Dirichlet_Esdo_exact_UV(T_cur0,rhs_cur0,...
                                Node,Node_flagx,Node_flagy,nnode,Discoef,uv_sol);   
else           % ... exact solution NOT known
    [T_cur,K] = Dirichlet_Esdo_matrix_UV(T_cur0,K0,Node_flagx,Node_flagy,nnode);
     rhs_cur      = Dirichlet_Esdo_rhs_UV(rhs_cur0,Node_flagx,Node_flagy,nnode);
end

%% Compute the purely elastic responce U_1 = U(t0), t0 = 0
uv_elast = T_cur\rhs_cur;   %  elastic step
uv_cur = uv_elast;

%% --------------------------------    Plot the current solution
UVu(1:nnode, 1) = uv_elast(1:nnode,1);
UVv(1:nnode, 1) = uv_elast(nnode+1:nallD,1);

% Plot the elastic solution 
switch test_problem
    case 0
%     nju = Discoef(1,Face_flag(1,1));
%     E   = Discoef(2,Face_flag(1,1))*Face_thick(1,1);
%     mju = E/(2*(1+nju));        lan = 2*mju*nju/(1-2*nju);
%     time_cur=0;
 %    [Uex,Vex] = UVtime_sol(Node(1,:)',Node(2,:)',time_cur);
      norm(UVu-Uex,inf),norm(UVv-Vex,inf)
      figure(4),clf,plot3(Node(1,:),Node(2,:),UVu(:,1)-Uex,'p')
      figure(5),clf,plot3(Node(1,:),Node(2,:),UVv(:,1)-Vex,'d')      
      % Test correctness of rhs
      figure(8),clf,plot(T_cur*uv_sol-rhs_cur,'r')
    case {9,10}
           figure(2),clf
   plot(Node(1,Surface_Nodes)*L_char, L_char*(Node(2,Surface_Nodes)+U_char*UVv(Surface_Nodes,1)'),'.')
   hold
   figure(3),clf
   plot(Node(1,Surface_Nodes)*L_char, L_char*(Node(2,Surface_Nodes)+U_char*UVu(Surface_Nodes,1)'),'.')
   hold
   disp(['Max displacement in x: ' num2str(max(abs(UVu(:,1)*U_char*L_char)))])
   disp(['Max displacement in y: ' num2str(max(abs(UVv(:,1)*U_char*L_char)))])
     figure(4),clf,plot3(Node(1,:),Node(2,:),UVu(:,1),'p')
     figure(5),clf,plot3(Node(1,:),Node(2,:),UVv(:,1),'d')
end
% --------------------------- End elastic step

% return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Viscoelastic solver - ODE DAE
%
% Assumptions:
% The positions of the Dirichlet b.c. do not change but the value of the
% values of the solution at the Dirichlet b.c. change with time.
%
% Forward Euler or Backward Euler for the time discretization
% Forward Euler: (currently not used)
%(K+A)ut_{k+1} = ((1-2*delta_t*alpha)*K-(1-delta_t*aplha)*A)ut_{k} + delta_t*alpha*f(k)
% Backward Euler : (currently used)
% ((1+2*alpha*delta_t)*K+(1+alpha*delta_t)A)U(k+1) = (K+A)U(K) + alpha*delta_t*f(k+1)
% Constant time step
    time = 0;
    flag_incr = 0; % ==1, incremental form

    delta_t_cur = min(delta_t_char,h);  % control the timestep

% I1: Form the matrix appearing on the right hane side (no b.c. applied)
tau = Maxwell_time_inv*delta_t_cur;
cK = 1+2*tau;
cA = 1+tau;
T_cur0 = cK*K0 + cA*A0;
S0         = K0 + A0;

% ----> Clear unnecessary arrays
clear A0

uv_til_cur = zeros(nall,1);  % initial condition for the auxiliary variable
% ------------------------------- Loop over time
% Implicit Euler
% Compute the load in the next time instance

cntr_cy = 1;
nexy = 0;
t_pr = [1,2,3,4,5,6,7,8,9,10,50,100];%,200,300,400,500,600,1000,[1500:500:10000]];
k = 1;
nface = size(Face_eorder9, 2);
UV=zeros(2*nnode,length(t_pr)+1);   % array to save the displacements in chosen years
next_store = 1;
UV(:,next_store) = uv_elast;                                 % store the elastic response

while (time_cur<=Tmax)
% %     disp(['Proceed with step ' int2str(k)])
% Set k=k+1, determine timestep dt(k) = t(k)-t(k-1)
%     delta_t_cur = delta_t_prev;  % take the same time step

% Save the previous solution
    uv_til_prev    = uv_til_cur;  % auxiliary varisble
    uv_prev          = uv_cur;      % displacements
 
 % Update the time
    time_cur = time_cur + delta_t_cur;

% Compute body forces at next time instance 
    rhs_b = spalloc(2*nnode,1,0); % no body forces
%     rhs_b = Assembly_ElAd_quadr_Q2_rhs(time_cur,Node,...
%                            Face_eorder9,Face_eorder4,...
%                            Face_flag,Face_thick,Discoef,...
%                            Gauss_point,Gauss_weight,...
%                            FUND_all,DERP_all,...
%                            nnode_lvl,nface_lvl,lvl_total,vec_coeff,wh);
 
% Compute surface forces  at next time instance 
       [rhs_s_new] = Assembly_rhs_Q2(Node,Edge,...
                               wh,l_ice,h_ice,rho_ice,grav,...
                               Load_Edges,Load_Edges_list,...
                               nnodeP,time_cur,T_BEG,T_LGM,T_EOG);
                           
    if flag_incr == 1
        rhs_s = rhs_s_new-rhs_s_prev;  % incremental surface force
    else
       rhs_s = rhs_s_new;                        % full furface force at next time step
    end
     % Update the rhs of the ODE
     rhs_cur0 = S0*uv_til_prev + tau*rhs_s;
 
     %% ------------------ impose Dirichlet b.c.; standard formulation ------------
if wh=='g0'    % ... exact solution known, homogeneous material assumed
    [Uex,Vex] = UVtime_sol(Node(1,:)',Node(2,:)',time_cur); 
    uv_sol = [Uex;Vex];
    [T_cur,rhs_cur] = Dirichlet_Esdo_exact_UV(T_cur0,rhs_cur0,...
                                Node,Node_flagx,Node_flagy,nnode,Discoef,uv_sol);   
else           % ... exact solution NOT known
    [T_cur,S] = Dirichlet_Esdo_matrix_UV(T_cur0,S0,Node_flagx,Node_flagy,nnode);
     rhs_cur      = Dirichlet_Esdo_rhs_UV(rhs_cur0,Node_flagx,Node_flagy,nnode);
end
 
    % Solve the system
    uv_til_cur = T_cur\rhs_cur;  
    
    % Recover uv_cur from uv_til_cur
    wrk = rhs_cur0-K0*uv_til_cur;
    wrk = Dirichlet_Esdo_rhs_UV(wrk,Node_flagx,Node_flagy,nnode);
    uv_cur = S\wrk;     
    
    if flag_incr == 1 % incremental form
        uv_cur = uv_elast + up_cur;
    end
    % Displacements
    UVu(1:nnode)  = uv_cur(1:nnode,1);
    UVv(1:nnode)  = uv_cur(nnode+1:2*nnode,1);

switch test_problem
    case {0,8}
    [Uex,Vex] = UVtime_sol(Node(1,:)',Node(2,:)',time_cur,nnodeP);
    figure(4),clf,plot3(Node(1,:)',Node(2,:)',UVu-Uex,'p')
    figure(5),clf,plot3(Node(1,:)',Node(2,:)',UVv-Vex,'d')

    norm(UVu-Uex),norm(UVv-Vex)
%     figure(4),clf,plot3(Node(1,:),Node(2,:),UVu(:,k),'p'),hold
%     figure(5),clf,plot3(Node(1,:),Node(2,:),UVv(:,k),'d'),hold
%
%     figure(4),plot3(Node(1,:),Node(2,:),UVu(:,k-1),'pc')
%     figure(5),plot3(Node(1,:),Node(2,:),UVv(:,k-1),'dc')
%
%     figure(4),plot3(Node(1,:),Node(2,:),Uex,'pr')
%     figure(5),plot3(Node(1,:),Node(2,:),Vex,'dr')
    case 10
      cy = floor(time_cur*T_char/secs_per_year);
      nexy = nexy + 1;
      [vmu,pmu]=max(abs(UVu(:)));
      [vmv,pmv]=max(abs(UVv(:)));
      save_maxu(1:2,nexy) = [cy;UVu(pmu)*U_char*L_char];
      save_maxv(1:2,nexy) = [cy;UVv(pmv)*U_char*L_char];
      if nexy>1
          time_dist = save_maxv(1,nexy)-save_maxv(1,nexy-1); % should be equal to 1 in this case
          if time_dist==0, time_dist = 1; end
%           if time_dist==0, keyboard,end
          sink_dist = save_maxv(2,nexy)-save_maxv(2,nexy-1);
           max_diff(nexy-1) = sink_dist/time_dist;
           disp(['Year ' num2str(cy) ', Average sinking per year: ' num2str(max_diff(nexy-1))])
      end

      if cy>=t_pr(cntr_cy)
          disp(['Year ' int2str(cy) ', Max displ in y and x: ' num2str(UVv(pmv)*U_char*L_char)  ' & ' num2str(UVu(pmu)*U_char*L_char)])
          next_store = next_store + 1;
          UV(:,next_store) = uv_cur;  % store the solution
          
          ydis = L_char*(U_char*reshape(UVv(id(:)), size(id)));
          xdis = L_char*(U_char*reshape(UVu(id(:)), size(id)));

          h2 = figure(2);
          p = plot(Node(1,id(1,:))*L_char, L_char*Node(2,id(1,:)) + ydis(1,:), '-.');
          h3 = figure(3);
          p = plot(Node(1,id(1,:))*L_char, L_char*Node(2,id(1,:)) + xdis(1,:), '-.');

        %   h4 = figure(4);
        %   imagesc(L_char*X(1,:), L_char*Y(:,1), ydis)
        %   set(gca,'YDir','normal')
        %   colormap(flipud(parula))
        %   colorbar
        %   title('Vertical displacement color map')
        %   ylabel('Depth')
        %   xlabel('Width')
          %
        %   h5 = figure(5);
        %   imagesc(L_char*X(1,:), L_char*Y(:,1), xdis)
        %   set(gca,'YDir','normal')
        %   colormap(flipud(parula))
        %   colorbar
        %   title('Horizontal displacement color map')
        %   ylabel('Depth')
        %   xlabel('Width')
          %
        %   h8 = figure(8);
        %   contour(L_char*X,L_char*Y,sigma_shear, linspace(min(sigma_shear(:)), max(sigma_shear(:)),50))
        %   colorbar
        %   title('Shear Stress contour')
        %   ylabel('Depth')
        %   xlabel('Width')

          drawnow
%           savefig(h2, './out/VerDispSurf')
%           savefig(h3, './out/HorDispSurf')
        %   savefig(h4, ['./out/VerDisp-y' num2str(cy)])
        %   savefig(h5, ['./out/HorDisp-y' num2str(cy)])
        %   savefig(h8, ['./out/ShearStress-y' num2str(cy)])

%           save(['./out/data-y' num2str(cy)], 'cy', 'L_char', 'U_char', 'M_char', 'Stress_cur', 'sigma_shear', 'Node', 'UVu', 'UVv', 'id', 'X', 'Y')
% %         save(['./out/data-y' num2str(cy)], 'cy', 'Pos', 'L_char', 'U_char', 'M_char', 'Stress_cur', 'Node', 'UVu', 'UVv', 'UVp', 'nnode_lvl', 'id', 'X', 'Y')

          cntr_cy = cntr_cy+1;
      end
end
     % Update norms
    norm_dif = norm(uv_cur-uv_prev);
    norm_cur = norm(uv_cur);
    norm_prev= norm(uv_prev);
%     disp(['norm(uv_cur-uv_prev) ' num2str(norm(uv_cur-uv_prev)) ...
%           ' ,norm(uv_cur) ' num2str(norm(uv_cur))])
end   % End Time loop
% max_diff = zeros(length(save_maxv),1);
% for nn=2:length(save_maxv)
%    max_diff(nn-1) = (save_maxv(2,nn)-save_maxv(2,nn-1))/(save_maxv(1,nn)-save_maxv(1,nn-1));
% end

return
