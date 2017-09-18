%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The viscoelastic solver
% ----- Initialize:
% % I1. Assembly the stiffness matrix,  S = A_elast + A_advec,
%     elasticity material coefficients at time tau=t
% I2. Assembly the elasticity stiffness matrix A_visc(0)
%     viscoelast coeff at time tau=t
% I3. Set k=1 add choose time t(k) and timestep dt(k)=t(k)-t(k-1)
% I4. Compute the initial rhs0
% I5. Choose initial guess for displacements and pressure, U(0)
% I6. Assembly the elasticity stiffness matrix A_visc(k)
%     viscoelast. coeff at time tau=dt(k)
% I7. W(k) = dt(k)/2*A_visc(k)*U(0)
%     rhs = rhs0 + W(0)
% I8. Solve (S - dk/2*A_visc(k))*U(k) = rhs;
% I9. Compute the integral [0->t(k)]:
%     I(K) = dk/2*A_visc(0)*U(k) + W(k)
% ----- Loop over time
% L1. Set k=k+1, determine timestep dt(k) = t(k)-t(k-1)
% L2. Assembly the elasticity stiffness matrix A_visc(k)
%     viscoelast coeff at time tau=tau+dk
% L3: Compute the rhs corresponding to a new load - rhs(k)
% L4. W(k) = dt(k)/2*A_visc(k)*U(k-1)
%     rhs  = rhs(k) + W(k) + exp(-dk/Maxwell_time)*I(k)
% L5. Solve (S - dt(k)/2*A_visc0)*U(k) = rhs;
% L6. Compute the integral term:
%     I(k)=I(k-1) + 1/2*dk(k-1)*A_visc(k)*U(k) + W(k)
% L7. Goto L1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I1. Assembly the global stiffness matrix
%     elasticity material coefficients at time tau=t
% K - elasticity part of S11
% A - advection part of S11  (S11=K+A)
%     [K+A     [B1 B2]^T ]
% S = [[B1 B2]    C      ]
% K and A - in separate displacement ordering
% C0 - pressure mass matix

%ilu_setup.droptol = actionILU(kkk);

% Precompute the basis functions for displacements and pressure
% and their derivatives (for the reference element, in the Gauss points)

%[Gauss_point,Gauss_weight] = Integr_weights_quad;
%Ga=[Gauss_point(1,:);Gauss_point(2,:);
%    Gauss_point(1,:).*Gauss_point(2,:);[1 1 1 1]];
%nip = size(Gauss_point,2);   % number of integration points
%FUND_all=cell(nip,1);
%DERD_all=cell(nip,1);
%for k=1:nip
%    FUND_all{k} = shape_fun_quad(Gauss_point,k);  % FUN(1,np)
%    DERD_all{k} = shape_der_quad(Gauss_point,k);  % DER(2,np),bilinear b.f.
%end

%%
%[Gauss_pointP,Gauss_weightP] = Integr_weights_quad;
[Gauss_point,Gauss_weight,nipD] = Integr_weights_quad_9;

% get the order of the nodes
[X,Y, id] = reorder_nodes(Node);

FUNP_all=cell(4,1);
DERP_all=cell(4,1);
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
nall  = nallD + nallP;

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
    %  vec_coeff = [0,-rho_earth*grav,0,0]
else
    disp('WITH PRE-STRESS ADVECTION.')  % with advection
    vec_coeff = [0,0,0,0]
    %   vec_coeff = [0,-rho_earth*grav,0, -rho_earth*grav]
    %  vec_coeff = [0,-rho_earth*grav,0,0]
end
% vec_coeff = vec_coeff.*advec_const  % scaling of the advection terms

theta = 0.5;    % in this case we use Trapetz method

k = 1;
secs_per_year = 365*24*3600; % Duration of a year in seconds 3.1536e+7
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[K0,A0,B1,B2,C0,...
    Face_estiffS,...
    rhs_db,rhs_pb] = Assembly_ElAd_quadrTH_Q2Q1(Node,Face_Node,Edge_Node',...
    Node_flagx,Node_flagy,...
    Face_Parent,Face_estiffS,...
    Face_eorder9,Face_eorder4,...
    Face_flag,Face_thick,Disco,Discoef,...
    Gauss_point,Gauss_weight,...
    FUNP_all,DERP_all,FUND_all,DERD_all,...
    nnode_lvl,nface_lvl,lvl_total,vec_coeff,wh); %pure elastic
rhs_b = [rhs_db;rhs_pb];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B0  = [B1;B2];
S_cur0 = [K0+A0 B0; B0' -C0];  % y-cordinate axis pointing upwards
Q_cur0 = [K0    B0; B0' -C0];
T_cur0 = S_cur0;
if (verbose ~= 0)
    disp('PRE-STRESS ADVECTION MATRIX WITH ''+''.')  %
end

%%% VERY UGLY FIX
global E nju mju lan rho
nju = Discoef(1,Face_flag(1,1));
E   = Discoef(2,Face_flag(1,1))*Face_thick(1,1);
mju = E/(2*(1+nju));        lan = 2*mju*nju/(1-2*nju);
rho = (1-2*nju)/(2*nju);    % mju/lan;


time_cur = 0;                      %  rhs & b.c. for the elastic response

%   Add surface forces
[rhs_ds,rhs_ps] = Assembly_rhsTH_Q2Q1(Node,Edge,...
    wh,l_ice,h_ice,rho_ice,grav,...
    Load_Edges,Load_Edges_list,...
    nnodeP,time_cur,T_BEG,T_LGM,T_EOG);

rhs_s = [rhs_ds;rhs_ps];
%  Compute the initial rhs0 (body forces + surface tractions)
rhs_cur0 = rhs_b + rhs_s;

%% ------------------------------------- impose Dirichlet b.c. ------------
if wh=='g0'    % ... exact solution known, homogeneous material assumed
    [Uex,Vex,Pex] = UVPtime_sol(Node(1,:)',Node(2,:)',time_cur,nnodeP);
    uvp_sol = [Uex;Vex;Pex];
    [T_cur,K,A,B,rhs_d,rhs_p] = Dirichlet_Esdo_exact(T_cur0,K0,A0,B0,rhs_cur0,...
        Node,Node_flagx,Node_flagy,nnode,nnodeP,Discoef,uvp_sol);
    rhs_cur=[rhs_d;rhs_p];
    
else           % ... exact solution NOT known
    [T_cur,K,A,B] = Dirichlet_Esdo_matrix(T_cur0,K0,A0,B0,Node_flagx,Node_flagy,nnode);
    [rhs_cur]   = Dirichlet_Esdo_rhs(rhs_cur0,Node_flagx,Node_flagy,nnode,nnodeP);
end
% Dirichlet_visco            %  impose b.c. both to S_cur and rhs_cur
% Q_cur = [K    B; B' -C0];
%% Compute the purely elastic responce U_1 = U(t0), t0 = 0
uvp_cur = T_cur\rhs_cur;   %  elastic step
%% --------------------------------    Plot the current solution
UVPu(1:nnode, 1) = uvp_cur(1:nnode,1);
UVPv(1:nnode, 1) = uvp_cur(nnode+1:nallD,1);
UVPp(1:nnodeP,1) = uvp_cur(nallD+1:nall,1);

figure(2),clf
plot(Node(1,Surface_Nodes)*L_char, (L_char*Node(2,Surface_Nodes)+L_char*U_char*UVPv(Surface_Nodes,1)'),'.')

hold
figure(3),clf
plot(Node(1,Surface_Nodes)*L_char, (L_char*Node(2,Surface_Nodes)+L_char*U_char*UVPu(Surface_Nodes,1)'),'.')
hold
disp(['Max displacenet in x: ' num2str(max(abs(UVPu(:,1)*U_char*L_char)))])
disp(['Max displacenet in y: ' num2str(max(abs(UVPv(:,1)*U_char*L_char)))])

% End elastic step
time = 0;
flag_incr = 1; % ==1, incremental form

if wh == 'g0', %(test_problem == 0)|(test_problem==9),
    nju = Discoef(1,Face_flag(1,1));
    E   = Discoef(2,Face_flag(1,1))*Face_thick(1,1);
    mju = E/(2*(1+nju));        lan = 2*mju*nju/(1-2*nju);
    rho = (1-2*nju)/(2*nju);    % mju/lan;
    time_cur=0;
    [Uex,Vex,Pex] = UVPtime_sol(Node(1,:)',Node(2,:)',time_cur,nnodeP);
    norm(UVPu-Uex),norm(UVPv-Vex), norm(UVPp-Pex)
    %     figure(4),clf,plot3(Node(1,:),Node(2,:),UVPu(:,1)-Uex,'p')
    %     figure(5),clf,plot3(Node(1,:),Node(2,:),UVPv(:,1)-Vex,'d')
    %     figure(6),clf,plot3(Node(1,1:nnodeP),Node(2,1:nnodeP),UVPp(:,1)-Pex,'o')
    %     figure(8),clf,plot(T_cur*uvp_sol),hold,plot(rhs_cur,'r')
    %     % Test correctness of rhs
    %     rhs_sol = T_cur*uvp_sol;
    %     www=rhs_sol-rhs_cur;
    %     [w1,w2]=find(www>1e-6);
    %     [rhs_sol(w1) rhs_d(w1) rhs_s(w1) rhs_cur(w1)]
end
%% ---  Checks for the solution in the case of uniform load
t_cur = 0;
bottom = min(Node(2,:)); %ymin
top        = max(Node(2,:)); %ymax
surface_load = rho_ice*grav*h_ice*rhs_const;
dnm = 1/(2*mju+lan);
dnms= 1/(mju*(2*mju/lan+1));
if vec_coeff(2)==0,  % no prestress
    exf = exp(Maxwell_time_inv*time_cur);
    gradV=-surface_load*dnms*exf;   % scaled pressure no prestress
    [max(UVPp) gradV max(UVPp)-gradV]
    UVPv_exact = surface_load*dnm*(bottom-Node(2,:))'*exf;
    UVPp_exact = gradV*ones(nnodeP,1);
else % ------------------- with prestress
    UVPv_exact = -rho_ice*h_ice/rho_earth*(exp(rho_earth*grav/(2*mju+lan)*(Node(2,:)-top)) - exp(rho_earth*grav/(2*mju+lan)*(bottom-top)))';
    UVPp_exact =  -surface_load*dnms*exp(rho_earth*grav/(2*mju+lan)*(Node(2,1:nnodeP)-top))';
end
figure(5),clf,plot3(Node(1,:),Node(2,:),UVPv_exact-U_char*UVPv,'v');
figure(6),clf,plot3(Node(1,1:nnodeP),Node(2,1:nnodeP),UVPp_exact-UVPp,'v');
% v=find(Node(2,:)==0); vp = find(Node(2,1:nnodeP)==0);
% sparse(UVPv(v)),sparse(UVPv_exact(v))
% [UVPp(vp), UVPp_exact(vp)]

fprintf('%9.7f & %1e & %e & %e \\\\ \n',nju,-max(abs(UVPv_exact))*L_char,-max(abs(UVPv))*L_char,norm((UVPv_exact-UVPv)*L_char,inf))
disp('End elastic step.')
return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Viscoelastic solver with accumulation of the visco-effects, Trapetz method

% Assumptions:
% The positions of the Dirichlet b.c. do not change but the value of the
% values of the solution at the Dirichlet b.c. change with time.
%
% Trapezoidal formula for numerical integration
% Constant time step
delta_t_cur = min(delta_t_char,h);  % control the timestep
% I1: Form the system matrix (no b.c. applied)
T_cur0 = S_cur0 - 0.5*Maxwell_time_inv*delta_t_cur*Q_cur0;

% ------------------------------- Loop over time
% uvp_sol_prev = [Uex;Vex;Pex];
cntr_cy = 1;
nexy = 0;
t_pr = [1,2,3,4,5,6,7,8,9,10,50,100];%,200,300,400,500,600,1000,[1500:500:10000]];
k = 1;
% Save the previous solution
uvp_prev    = uvp_cur;
uvp_elast   = uvp_cur;
rhs_s_prev  = rhs_s;
norm_diff = 1000;
disp(['Norm rhs: ' num2str(norm(rhs_cur))])

nface = size(Face_eorder9, 2);
Stress_cur    = zeros(nface, 3);
Stress_e_cur  = zeros(nface, 3);
Stress_prev   = zeros(nface, 3);
Stress_e_prev = zeros(nface, 3);
while (time_cur<=Tmax) || (norm_diff>1e-5)
    % %     disp(['Proceed with step ' int2str(k)])
    %     uvp_sol_prev= uvp_sol;
    % Set k=k+1, determine timestep dt(k) = t(k)-t(k-1)
    % %     delta_t_cur = delta_t_prev;  % take the same time step
    
    Stress_prev  = Stress_cur;
    Stress_e_prev = Stress_e_cur;
    
    % Recover stress field from displacements
    % Stress_cur = SStime_app(UVPu,UVPv,Node,Edge,Face,Face_eorder9,Face_Parent,...
    %  nnode,hx,hy,time_cur,Maxwell_time_inv);
    % [Pos, Stress_cur] = SStime_midpoint(UVPu,UVPv,UVPp, Node,Edge,Face,Face_eorder9,Face_Parent,...
    %                                              nnode,hx,hy,time_cur,Maxwell_time_inv);
    
    if k == 1
        [Pos, Stress_cur] = SStime_midpoint(UVPu,UVPv,UVPp, Node,Edge,Face,Face_eorder9,Face_Parent,...
            nnode,hx,hy,time_cur,Maxwell_time_inv);
    else
        [Stress_cur, Stress_e_cur] = SStime_mid_iter(Stress_prev, Stress_e_prev, delta_t_cur, UVPu,UVPv,UVPp, Node,Edge,Face,Face_eorder9,Face_Parent,...
            nnode,hx,hy,time_cur,Maxwell_time_inv,k);
    end
    
    
    % Update the time
    
    time_cur = time_cur + delta_t_cur;
    k = k + 1;
    
    % Compute the rhs corresponding to a new load - rhs(k)
    % Compute the body forces
    [rhs_db,rhs_pb] = Assembly_ElAd_quadrTH_Q2Q1_rhs(time_cur,Node,...
        Face_eorder9,Face_eorder4,...
        Face_flag,Face_thick,Discoef,...
        Gauss_point,Gauss_weight,...
        FUND_all,DERP_all,...
        nnode_lvl,nface_lvl,lvl_total,vec_coeff,wh);
    rhs_b = [rhs_db; rhs_pb];
    % Compute the surface forces
    if flag_incr == 1
        rhs_s = zeros(nall,1);
    else
        [rhs_ds,rhs_ps] = Assembly_rhsTH_Q2Q1_cur(Node,Edge,Stress_cur,...
            wh,l_ice,h_ice,rho_ice,grav,...
            Load_Edges,Load_Edges_list,...
            Top_Edges,Top_Edges_list,...
            Right_Edges,Right_Edges_list,...
            nnodeP,time_cur,T_BEG,T_LGM,T_EOG);
        rhs_s = [rhs_ds;rhs_ps];
    end
    % Update the memory term  and the visco boundary term
    if k == 2
        W  = 0.5*Maxwell_time_inv*delta_t_cur*exp(-Maxwell_time_inv*delta_t_cur)*Q_cur0*uvp_prev;
        WV = 0.5*Maxwell_time_inv*delta_t_cur*(rhs_s+exp(-Maxwell_time_inv*delta_t_cur)* rhs_s_prev);
    else
        W  =  exp(-Maxwell_time_inv*delta_t_cur)*W  +  ...
            0.5*Maxwell_time_inv*delta_t_cur*exp(-Maxwell_time_inv*delta_t_cur)*Q_cur0*uvp_prev ;
        WV =  exp(-Maxwell_time_inv*delta_t_cur)*WV  +  ...
            0*0.5*Maxwell_time_inv*delta_t_cur*(rhs_s+exp(-Maxwell_time_inv*delta_t_cur)*rhs_s_prev);
    end
    % Save the previous solution
    uvp_prev    = uvp_cur;
    rhs_s_prev  = rhs_s;
    % Add the memory term
    rhs_cur0    = rhs_b + rhs_s + W - WV;  %
    disp(['Norm rhs_s:                        '  num2str(norm(rhs_s))])
    disp(['Norm rhs_W:                        '  num2str(norm(W))])
    disp(['Norm rhs_WV:                       '  num2str(norm(WV))])
    disp(['Norm rhs_cur before b.c.: '  num2str(norm(rhs_cur0))])
    
    % Apply Dirichlet b.c., corresponding to time_cur
    if wh=='g0'    % ... exact solution known, homogeneous material assumed
        [Uex,Vex,Pex] = UVPtime_sol(Node(1,:)',Node(2,:)',time_cur,nnodeP);
        uvp_sol = [Uex;Vex;Pex];
        [T_cur,K,A,B,rhs_d,rhs_p] = Dirichlet_Esdo_exact(T_cur0,K0,A0,B0,rhs_cur0,...
            Node,Node_flagx,Node_flagy,nnode,nnodeP,Discoef,uvp_sol);
        rhs_cur = [rhs_d;rhs_p];
        
    else           % ... exact solution NOT known
        [T_cur,K,A,B] = Dirichlet_Esdo_matrix(T_cur0,K0,A0,B0,Node_flagx,Node_flagy,nnode);
        rhs_prev      = rhs_cur;
        [rhs_cur]     = Dirichlet_Esdo_rhs(rhs_cur0,Node_flagx,Node_flagy,nnode,nnodeP);
        %        disp(['Norm rhs_cur after  b.c.: '  num2str(norm(rhs_cur))])
    end
    
    switch test_problem
        
        case {0,8}
            uvp_cur = T_cur\rhs_cur;
        case {9,10}
            if flag_incr == 1 % incremental form
                uvp_cur_incr = T_cur\rhs_cur;
                uvp_cur = uvp_elast + uvp_cur_incr;
            else
                uvp_cur = T_cur\rhs_cur;
            end
    end
    % Check the relation between div u and p
    %   figure(9),plot(uvp_cur(1:2*nnode,1)),drawnow
    %   figure(10),plot(uvp_cur(1+2*nnode:end,1)),drawnow,pause(2)
    %   wrk=T_cur(2*nnode+1:end,1:2*nnode)*uvp_cur(1:2*nnode,1) +T_cur(2*nnode+1:end,2*nnode+1:end)*uvp_cur(2*nnode+1:end);
    % Update the memory integral
    W = W + 0.5*Maxwell_time_inv*delta_t_cur*Q_cur0*uvp_cur ;
    
    % Displacements
    UVPu(1:nnode)  = uvp_cur(1:nnode,1);
    UVPv(1:nnode)  = uvp_cur(nnode+1:2*nnode,1);
    UVPp(1:nnodeP) = uvp_cur(2*nnode+1:2*nnode+nnodeP,1);
    
    switch test_problem
        case {0,8}
            [Uex,Vex,Pex] = UVPtime_sol(Node(1,:)',Node(2,:)',time_cur,nnodeP);
            figure(4),clf,plot3(Node(1,:)',Node(2,:)',UVPu-Uex,'p')
            figure(5),clf,plot3(Node(1,:)',Node(2,:)',UVPv-Vex,'d')
            figure(6),clf,plot3(Node(1,1:nnodeP)',Node(2,1:nnodeP)',UVPp-Pex,'o')
            
            norm(UVPu-Uex),norm(UVPv-Vex), norm(UVPp-Pex)
            %     figure(4),clf,plot3(Node(1,:),Node(2,:),UVPu(:,k),'p'),hold
            %     figure(5),clf,plot3(Node(1,:),Node(2,:),UVPv(:,k),'d'),hold
            %     figure(6),clf,plot3(Node(1,1:nnodeP),Node(2,1:nnodeP),UVPp(:,k),'o') ,hold
            %
            %     figure(4),plot3(Node(1,:),Node(2,:),UVPu(:,k-1),'pc')
            %     figure(5),plot3(Node(1,:),Node(2,:),UVPv(:,k-1),'dc')
            %     figure(6),plot3(Node(1,1:nnodeP),Node(2,1:nnodeP),UVPp(:,k-1),'oc')
            %
            %     figure(4),plot3(Node(1,:),Node(2,:),Uex,'pr')
            %     figure(5),plot3(Node(1,:),Node(2,:),Vex,'dr')
            %     figure(6),plot3(Node(1,1:nnodeP),Node(2,1:nnodeP),Pex,'or')
        case {9,10}
            cy = floor(time_cur*T_char/secs_per_year);
            nexy = nexy + 1;
            [vmu,pmu]=max(abs(UVPu(:)));
            [vmv,pmv]=max(abs(UVPv(:)));
            save_maxu(1:2,nexy) = [cy;UVPu(pmu)*U_char*L_char];
            save_maxv(1:2,nexy) = [cy;UVPv(pmv)*U_char*L_char];
            %       figure(4),clf,plot3(Node(1,:)',Node(2,:)',UVPu,'p')
            %       figure(5),clf,plot3(Node(1,:)',Node(2,:)',UVPv,'d')
            %       figure(6),clf,plot3(Node(1,1:nnodeP)',Node(2,1:nnodeP)',UVPp,'o')
            if nexy>1
                time_dist = save_maxv(1,nexy)-save_maxv(1,nexy-1); % should be equal to 1 in this case
                if time_dist==0, time_dist = 1; end
                %           if time_dist==0, keyboard,end
                sink_dist = save_maxv(2,nexy)-save_maxv(2,nexy-1);
                max_diff(nexy-1) = sink_dist/time_dist;
                disp(['Year ' num2str(cy) ', Average sinking per year: ' num2str(max_diff(nexy-1))])
            end
            
            if cy>=t_pr(cntr_cy)
                disp(['Year ' int2str(cy) ', Max displ in y and x: ' num2str(UVPv(pmv)*U_char*L_char)  ' & ' num2str(UVPu(pmu)*U_char*L_char)])
                %           UVPv(pmv)*U_char*L_char
                %           UVPu(pmu)*U_char*L_char,wait
                
                ydis = L_char*(U_char*reshape(UVPv(id(:)), size(id)));
                xdis = L_char*(U_char*reshape(UVPu(id(:)), size(id)));
                
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
                
                %           save(['./out/data-y' num2str(cy)], 'cy', 'L_char', 'U_char', 'M_char', 'Stress_cur', 'sigma_shear', 'Node', 'UVPu', 'UVPv', 'id', 'X', 'Y')
                save(['./out/data-y' num2str(cy)], 'cy', 'Pos', 'L_char', 'U_char', 'M_char', 'Stress_cur', 'Node', 'UVPu', 'UVPv', 'UVPp', 'nnode_lvl', 'id', 'X', 'Y')
                
                cntr_cy = cntr_cy+1;
            end
    end
    % Update norms
    norm_dif = norm(uvp_cur-uvp_prev);
    norm_cur = norm(uvp_cur);
    norm_prev= norm(uvp_prev);
    %     disp(['norm(uvp_cur-uvp_prev) ' num2str(norm(uvp_cur-uvp_prev)) ...
    %           ' ,norm(uvp_cur) ' num2str(norm(uvp_cur))])
    %% ---  Checks for the solution in the case of uniform load
    mjup = mju*exp(Maxwell_time_inv*t_cur);
    lanp = lan*exp(Maxwell_time_inv*t_cur);
    dnm = 1/(2*mjup+lanp);
    dnms= 1/(mjup*(2*mjup/lanp+1));
    if vec_coeff(2)==0,  % no prestress
        exf = exp(Maxwell_time_inv*time_cur);
        gradV=-surface_load*dnms*exf;   % scaled pressure no prestress
        [max(UVPp) gradV max(UVPp)-gradV]
        UVPv_exact = surface_load*dnm*(bottom-Node(2,:))'*exf;
        UVPp_exact = gradV*ones(nnodeP,1);
    else % ------------------- with prestress
        UVPv_exact = -rho_ice*h_ice/rho_earth*(exp(rho_earth*grav/(2*mjup+lanp)*(Node(2,:)-top)) - exp(rho_earth*grav/(2*mjup+lanp)*(bottom-top)))';
        UVPp_exact =  -surface_load*dnms*exp(rho_earth*grav/(2*mjup+lanp)*(Node(2,1:nnodeP)-top))';
    end
    figure(5),clf,plot3(Node(1,:),Node(2,:),UVPv_exact-UVPv,'v');
    figure(6),clf,plot3(Node(1,1:nnodeP),Node(2,1:nnodeP),UVPp_exact-UVPp,'v');
end   %
% max_diff = zeros(length(save_maxv),1);
% for nn=2:length(save_maxv)
%    max_diff(nn-1) = (save_maxv(2,nn)-save_maxv(2,nn-1))/(save_maxv(1,nn)-save_maxv(1,nn-1));
% end

return
