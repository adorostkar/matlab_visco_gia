% Given: Node,Edge,Face
% asembling of the time-dependent body force, per timestep

function [rhs_d,rhs_p] = Assembly_ElAd_quadrTH_Q2Q1_rhs(time,Node,...
        Face_eorder9,Face_eorder4,...
        Face_flag,Face_thick,Discoef,...
        Gauss_point,Gauss_weight,...
        FUND_all,DERP_all,...
        nnode_lvl,nface_lvl,lvl_total,vec_coeff,wh)
    global verbose
    global mju lan
    
    nnodeD = nnode_lvl(lvl_total);   %number of nodes for the displacements
    nnodeP = nnode_lvl(lvl_total-1); %number of nodes for the pressure
    nface  = nface_lvl(lvl_total-1); %number of faces for the pressure
    dim   = 2;	   %
    dofd =  9;
    dofD = dofd*dim; % ndof = nip*dim
    dofP =  4;
    
    nrmmax = 0;
    if(verbose ~= 0)
        disp('Begin allocating memory...')
    end
    
    rhs_d= zeros(2*nnodeD,1);
    rhs_p= zeros(nnodeP,1);
    nallV = 2*nnodeD;
    nallP = nnodeP;
    
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
        
        % - - - - - generation of the element matrices:
        
        M_elem=Assm_ElAdSaddle_quadrTH_Q2Q1_rhs(Gauss_point,Gauss_weight,...
            FUND_all,DERP_all,CoordP);
        
        macro_node_listD=[local_node_listD;local_node_listD+nnodeD]; % sdo
        ZZ = zeros(size(M_elem));
        if ~strcmp(wh,'g0')
            % % No body forces in this case
            rhs_elemD = zeros(dofD,1);
        else
            bforce = bodyf_time(Node(1,local_node_listD)',Node(2,local_node_listD)',time,...
                mju,lan,vec_coeff);
            
            % % - - - - - bforce is a (ndof x 1) vector (for quadrilaterals)
            % %           IN separate displacements ordering !!!
            MM = [M_elem ZZ;ZZ M_elem];
            rhs_elemD = MM*bforce;
        end
        
        % Assembly of the global rhs
        for ii=1:dofD,
            is = macro_node_listD(ii);
            rhs_d(is,1) = rhs_d(is,1) + rhs_elemD(ii,1);
        end
        
    end             % end pressure iface_p-loop
    return
    
    
