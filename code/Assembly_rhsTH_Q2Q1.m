% Given: Node,Edge
% asembling of the right-hand-side vector (surface load only)
% The load is time-dependent
%
%  o-----*-----o
%  |                |
%  |                |
%  *       *      *
%  |                |
%  |                |
%  o-----*-----o
%

function [rhs_sd,rhs_sp] = Assembly_rhsTH_Q2Q1(Node,Edge,...
        wh,l_ice,h_ice,rho_ice,grav,...
        Load_Edges,Load_Edges_list,...
        nnodeP,time_cur,T_BEG,T_LGM,T_EOG)
    global rhs_const h
    global test_problem Maxwell_time_inv
    global mju lan
    global rho_ice grav h_ice
    
    nnode = size(Node,2);   % number of nodes
    
    rhs_sd = zeros(2*nnode,1);
    rhs_sp = zeros(nnodeP,1);
    
    switch test_problem
        case  {0,8}
            % - - - - - - no surface forces to add
            return
        case {9,10}
            % - - - - - - Add surface forces to the right-hand side
            %Walk on all edges, where traction force is applied,
            %in this case, edges under the ice
            % For each such edge, integrate the three basis functions
            % Integrate numerically along the boundary edge (three point quadrature)
            
            % We use the Simpson quadrature formula
            % On the interval [-1,1] the integral of the three b.f. is 2
            
            for iiedge = 1:length(Load_Edges_list),
                iedge = Load_Edges_list(iiedge);
                cnodes= Load_Edges(:,iiedge);
                fnodes= Edge(:,iedge);
                % Find the midpoint, belonging to this (coarse) edge
                midnode=max(Edge(:,iedge));
                nn = [min(cnodes),max(cnodes),midnode];
                Coord = Node(:,Load_Edges(:,iiedge))';  % Coord(x/y,no_point)
                
                %      for k=1:2, nn(k)=find(local_node_list==Load_Edges(k,iiedge));end
                %            nn(3)=find(local_node_list==max(Edge(:,iedge)));
                %
                xmid = Node(1,midnode);
                ymid = Node(2,midnode);
                xdif = abs(diff(Coord(:,1)));
                ydif = abs(diff(Coord(:,2)));
                len  = sqrt(xdif^2+ydif^2);  % the length of the edge
                
                Int = len/6*[1;1;4];  % <------ Simpson's formula
                if test_problem >= 10, % the GIA problem
                    h_ice_cur = height_ice(wh,xmid, ymid, h_ice,time_cur,T_BEG,T_LGM,T_EOG);
                    n=[0,1]; % outer unit normal on the top line
                    G1 = zeros(3,1);
                    %             G2 = -grav*rho_ice*h_ice_cur*(2-exp(-Maxwell_time_inv*time_cur))*Int;
                    G2 = -grav*rho_ice*h_ice_cur*Int;
                    G2 = rhs_const*G2;
                else
                    G1 = zeros(3,1);
                    uniform_load = -rho_ice*grav*h_ice*Int; % unit square
                    %              h_ice_cur = height_ice(wh,xmid, ymid, h_ice,time_cur,T_BEG,T_LGM,T_EOG);
                    %              uniform_load = -grav*rho_ice*h_ice_cur*Int;
                    G2 = rhs_const*uniform_load;
                end
                % ---------------------------------------------
                rhs_sd(nn,1)       = rhs_sd(nn,1)       + G1;
                rhs_sd(nn+nnode,1) = rhs_sd(nn+nnode,1) + G2;
            end   % end loop over edges
    end   % end case
    
    return
