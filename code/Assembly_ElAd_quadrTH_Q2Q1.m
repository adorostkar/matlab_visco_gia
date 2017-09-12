% Given: Node,Edge,Face
% asembling of the global stiffness matrix (Neumann everywhere)
% K(nnode,nnode) and the right-hand-side vector (only body forces)
% K is in saddle-point form; kinematic pressure intruduced
% Stable discretization (Q2-Q1 Taylor-Hood elements) !!!
% --> C  is the 22 matrix block



function [K,A,B1,B2,C,Face_estiffS,rhs_d,rhs_p] = ...
              Assembly_ElAd_quadrTH_Q2Q1(Node,Face_Node,Node_Edge,...
                                         Node_flagx,Node_flagy,...
                                         Face_Parent,Face_estiffS,...
                                         Face_eorder9,Face_eorder4,...
                                         Face_flag,Face_thick,Disco,Discoef,...
                                         Gauss_point,Gauss_weight,...
                                         FUNP_all,DERP_all,FUND_all,DERD_all,...
                                         nnode_lvl,nface_lvl,lvl_total,vec_coeff,wh)
    global verbose
    global advec_const
    global lan mju

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
    B1 = spalloc(nnodeD,nnodeP,5*nnodeD); % (div velocity, pressure test func)
    B2 = spalloc(nnodeD,nnodeP,5*nnodeD);
    C  = spalloc(nnodeP,nnodeP,9*nnodeP); % mass - pressure
    S  = spalloc(2*nnodeD+nnodeP,2*nnodeD+nnodeP,9*(2*nnodeD+nnodeP)); % saddle
    rhs_d= zeros(2*nnodeD,1);
    rhs_p= zeros(nnodeP,1);
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
    nextK = 0;

    lengthA=nface*dofD*dofD;
    AI  = zeros(lengthA,1);
    AJ  = zeros(lengthA,1);
    AV  = zeros(lengthA,1);
    nextA = 0;

    lengthB1=nface*dofD*dofP;     % B1 -> mju*p*phi_x
    B1I = zeros(lengthB1,1);
    B1J = zeros(lengthB1,1);
    B1V = zeros(lengthB1,1);
    nextB1 = 0;

    lengthB2=nface*dofD*dofP;    % B2 -> mju*p*phi_y
    B2I = zeros(lengthB2,1);
    B2J = zeros(lengthB2,1);
    B2V = zeros(lengthB2,1);
    nextB2 = 0;

    lengthC=nface*dofP*dofP;     % the _22 block (mass matrix)
    CI  = zeros(lengthC,1);
    CJ  = zeros(lengthC,1);
    CV  = zeros(lengthC,1);
    nextC = 0;

    if(verbose ~= 0)
        disp('...end allocating memory.')
    end

    for iface_p=1:nface,     % ---> walk on the parent (pressure) elements
        %    disp(['Pressure face ' int2str(iface_p)])
        local_node_listP = Face_eorder4(:,iface_p);
        local_node_listD0= Face_eorder9(:,iface_p);
        local_node_listD=order_face_nodes_quadQ2(Node,local_node_listP,local_node_listD0);

        CoordP(1:4,1) = Node(1,local_node_listP)';  % Coord4(4,2)
        CoordP(1:4,2) = Node(2,local_node_listP)';
        CoordD(1:9,1) = Node(1,local_node_listD)';  % Coord9(9,2)
        CoordD(1:9,2) = Node(2,local_node_listD)';


        % % - - - - - generation of the pressure element matrices:
        %     [C_elem,M_elem0]=Assm_quadrTH(Gauss_weightP,...
        % 			               FUNP_all,DERP_all,CoordP,wh);
        %
        % ------------------------ determine nju, E
        nju = Discoef(1,Face_flag(iface_p,1));
        E   = Discoef(2,Face_flag(iface_p,1))*Face_thick(iface_p,1);

        % ------------------------ compute the coefficients lan, mju
        mju = E/(2*(1+nju));        lan = 2*mju*nju/(1-2*nju);
        rho = (1-2*nju)/(2*nju);    % mju/lan;

        % - - - - - generation of the element matrices:

        [El_elem,Ad_elem,B1_elem,B2_elem,M_elem,M0_elem, L_elem9,R_elem9]=...
            Assm_ElAdSaddle_quadrTH_Q2Q1(Gauss_point,Gauss_weight,...
                                         FUNP_all,DERP_all,FUND_all,DERD_all,...
                                         CoordP,CoordD,vec_coeff,nju,wh);

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
        if ~strcmp(wh,'g0')
            % % No body forces in this case
            rhs_elemD = zeros(dofD,1);
        else
            bforce = bodyf_time(Node(1,local_node_listD)',Node(2,local_node_listD)',...
                                0,mju,lan,vec_coeff);

            % % - - - - - bforce is a (ndof x 1) vector (for quadrilaterals)
            % %           IN separate displacements ordering !!!
            MM = [M_elem ZZ;ZZ M_elem];
            rhs_elemD = MM*bforce;
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
                KV(nextK) = mju*El_elem(ii,jj);

                LI(nextK) = is;
                LJ(nextK) = js;
                LV(nextK) = L_elem9(ii,jj);
                RI(nextK) = is;
                RJ(nextK) = js;
                RV(nextK) = R_elem9(ii,jj);

                nextA = nextA + 1;
                AI(nextA) = is;
                AJ(nextA) = js;
                AV(nextA) = Ad_elem(ii,jj); %all advection coeffs must already
                %be included in the element matrix
            end
            rhs_d(is,1) = rhs_d(is,1) + rhs_elemD(ii,1);
        end

        for ii=1:dofd,
            is = local_node_listD(ii);
            for jj=1:dofP,
                js = local_node_listP(jj);
                % 	  B1(is,js)= B1(is,js)+ mju*B1_elem(ii,jj);
                nextB1 = nextB1 + 1;
                B1I(nextB1) = is;
                B1J(nextB1) = js;
                B1V(nextB1) = mju*B1_elem(ii,jj);

                %	  B2(is,js)= B2(is,js)+ mju*B2_elem(ii,jj);
                nextB2 = nextB2 + 1;
                B2I(nextB2) = is;
                B2J(nextB2) = js;
                B2V(nextB2) = mju*B2_elem(ii,jj);
            end
        end

        for ii=1:dofP,
            is = local_node_listP(ii);
            for jj=1:dofP,
                js = local_node_listP(jj);

                %         C(is,js) = C(is,js) - rho*mju*M_elem0(ii,jj);
                nextC = nextC + 1;
                CI(nextC) = is;
                CJ(nextC) = js;
                CV(nextC) = rho*mju*M0_elem(ii,jj);
            end
        end
        %     B_elem = [B1_elem;B2_elem];
        %     S_macro = [2*mju*El_elem-advec_const*Ad_elem mju*B_elem; mju*B_elem' -rho*mju*M0_elem];
        %     Face_estiffS(1:22,1:22,iface_p) = S_macro;
        %
        %     Q11=2*mju*El_elem-advec_const*Ad_elem;
        %     Q11I=inv(Q11+0.05*eye(18));
        %     QQ=mju*B_elem'*Q11I*mju*B_elem;
        %     EM=eigs(M0_elem,QQ,1);
        %     if EM>nrmmax, nrmmax=EM;end

    end             % end pressure iface_p-loop

    KI = KI(1:nextK,1); KJ = KJ(1:nextK,1); KV = KV(1:nextK,1);
    AI = AI(1:nextA,1); AJ = AJ(1:nextA,1); AV = AV(1:nextA,1); %
    B1I=B1I(1:nextB1,1);B1J=B1J(1:nextB1,1);B1V=B1V(1:nextB1,1);
    B2I=B2I(1:nextB2,1);B2J=B2J(1:nextB2,1);B2V=B2V(1:nextB2,1);

    % W1I=W1I(1:nextW1,1);W1J=W1J(1:nextW1,1);W1V=W1V(1:nextW1,1); % B^T
    % W2I=W2I(1:nextW2,1);W2J=W2J(1:nextW2,1);W2V=W2V(1:nextW2,1);

    CI = CI(1:nextC,1); CJ = CJ(1:nextC,1); CV = CV(1:nextC,1);

    K  = sparse(KI, KJ, KV);
    A  = advec_const*sparse(AI, AJ, AV);
    B1 = sparse(B1I,B1J,B1V);
    B2 = sparse(B2I,B2J,B2V);
    C  = sparse(CI, CJ, CV);
    return

    % B = spalloc(2*nnodeD,nnodeP,0);
    % C = spalloc(nnodeP,nnodeP,0);
    % LI = LI(1:nextK,1); LJ = LJ(1:nextK,1); LV = LV(1:nextK,1); % Laplace only
    % RI = RI(1:nextK,1); RJ = RJ(1:nextK,1); RV = RV(1:nextK,1); % rot only
    % L  = 2*mju*sparse(LI, LJ, LV); % Laplace only
    % R  = mju*sparse(RI, RJ, RV); % rot only
    % save Laplace_mat L
    % save rot_mat R
    % save advec_mat A
    % save mass_p_mat C
    % save B1_mat B1
    % save B2_mat B2
    %
    % % K-(2*mju*L+mju*R)=0, tested
    %
    % % ------
