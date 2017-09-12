%%%%%%%%%%%%%%%%%%%% Easy DD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------- Uniform efinement of quadrilaterals --------------
% This routine refines the very coarse mesh only once
% to mark properly the nodes and edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Node,Edge,Face,...
        Node_flagx,Node_flagy,Edge_flagx,Edge_flagy,Face_flag,...
        Face_thick,...
        Face_Node,Face_Parent,Face_eorderD,Face_eorderP,...
        nface_lvl,nnode_lvl,lvl] = ...
            my_Refine_quadr_hier(Node,Edge,Face,...
                                 Node_flagx,Node_flagy,...
                                 Edge_flagx,Edge_flagy,Face_flag,Face_thick,...
                                 nface_lvl,nnode_lvl,lvl)

    nnode = size(Node,2);
    nedge = size(Edge,2);
    nface = size(Face,2);

    new_node = 0;
    new_edge = 0;

    Face_Parent   = zeros(4,nface,1);
    Face_Node     = zeros(4,4*nface,2);
    Face_eorderD = zeros(9,nface);
    Face_eorderP = zeros(4,nface);

    % Make Node, Edge and Face immediately as large as needed
    % Every face will be divided in 4
    Face(1,4*nface)      = 0;
    Face_flag(4*nface,1) = 0;
    Face_thick(4*nface,1)= 1;
    % every edge will produce a new node and every face - one internal node
    Node(1,nnode+nedge+nface)       = 0;
    Node_flagx(nnode+nedge+nface,1) = 0;
    Node_flagy(nnode+nedge+nface,1) = 0;
    % Euler #F + #N - #E = 2  but outside region is also a Face
    Edge(1,4*nface+2*nedge)       = 0;
    Edge_flagx(4*nface+2*nedge,1) = 0;
    Edge_flagy(4*nface+2*nedge,1) = 0;

    % -------------------------- halve the edges
    for iedge = 1:nedge,
        inode1 = Edge(1,iedge); inode2 = Edge(2,iedge);
        %    disp(['Edge ',int2str(iedge)])
        %disp(['End_1 ' int2str(inode1) ', coord: [' num2str(x1) ',' num2str(y1) ']'])
        %disp(['End_2 ' int2str(inode1) ', coord: [' num2str(x2) ',' num2str(y2) ']'])
        %disp(['Midpoint coord: [' num2str(xM) ',' num2str(yM) ']']),wait
        new_node = new_node + 1;
        Node(:,nnode+new_node) = (Node(:,inode1)+Node(:,inode2))*0.5;
        new_edge = new_edge + 1;
        Edge(1,nedge+new_edge) = nnode+new_node;
        Edge(2,nedge+new_edge) = Edge(2,iedge);
        Edge(2,iedge) = nnode+new_node;
        Edge_flagx(nedge+new_edge,1) = Edge_flagx(iedge,1);
        Edge_flagy(nedge+new_edge,1) = Edge_flagy(iedge,1);
        Node_flagx(nnode+new_node,1) = Edge_flagx(iedge,1);
        Node_flagy(nnode+new_node,1) = Edge_flagy(iedge,1);
    end
    % --------------------- add internal edges and new faces

    new_face  = nface+1;
    iedge     = zeros(1,6);
    node_mid  = zeros(1,3);
    Face2Edge = zeros(1,3);
    CurEdge   = nedge*2+1;   % CurEdge=size(Edge,2)+1;

    % ---------------- The ONLY assumption made is that the order of the edge
    %                  numbers defining a face is such that they make a
    %                  closed path (kept through the refinement)

    for iface = 1:nface,  % ------------- begin loop over edges
        for j = 1:4,
            old_edge(j) = Face(j,iface);
            new_edge(j) = old_edge(j) + nedge; % similar as Edge_temp
        end
        % ----------------------- make an oriented local copy of the edges and nodes
        local_edge = zeros(8,1);
        local_Edge = zeros(2,8);
        local_edge(1)   = old_edge(1);
        local_Edge(:,1) = Edge(:,old_edge(1));
        local_edge(2)   = new_edge(1);
        local_Edge(:,2) = my_flip(Edge(:,new_edge(1)),local_Edge(2,1),1);
        n = 2;
        for j=2:4,
            if ( local_Edge(2,n)==Edge(1,old_edge(j)) ) | ...
                    ( local_Edge(2,n)==Edge(2,old_edge(j)) ),
                % if the last point belongs to the "old" part of the next edge...
                local_edge(n+1)   = old_edge(j);
                local_Edge(:,n+1) = my_flip(Edge(:,old_edge(j)),local_Edge(2,n),1);
                local_edge(n+2)   = new_edge(j);
                local_Edge(:,n+2) = my_flip(Edge(:,new_edge(j)),local_Edge(2,n+1),1);
            elseif ( local_Edge(2,n)==Edge(1,new_edge(j)) ) | ...
                    ( local_Edge(2,n)==Edge(2,new_edge(j)) ),
                % if the last point belongs to the "new" part of the next edge...
                local_edge(n+1)   = new_edge(j);
                local_Edge(:,n+1) = my_flip(Edge(:,new_edge(j)),local_Edge(2,n),1);
                local_edge(n+2)   = old_edge(j);
                local_Edge(:,n+2) = my_flip(Edge(:,old_edge(j)),local_Edge(2,n+1),1);
            else
                local_edge(3:n*2)  =local_edge(1:(n-1)*2);
                local_Edge(:,3:n*2)=local_Edge(:,1:(n-1)*2);
                if ( local_Edge(1,3)==Edge(1,old_edge(j)) ) | ...
                        ( local_Edge(1,3)==Edge(2,old_edge(j)) ),
                    % if the first point belongs to the "old" part of the next edge...
                    local_edge(2)   = old_edge(j);
                    local_Edge(:,2) = my_flip(Edge(:,old_edge(j)),local_Edge(1,3),2);
                    local_edge(1)   = new_edge(j);
                    local_Edge(:,1) = my_flip(Edge(:,new_edge(j)),local_Edge(1,2),2);
                elseif   ( local_Edge(1,3)==Edge(1,new_edge(j)) ) | ...
                        ( local_Edge(1,3)==Edge(2,new_edge(j)) ),
                    % if the first point belongs to the "new" part of the next edge...
                    local_edge(2)   = new_edge(j);
                    local_Edge(:,2) = my_flip(Edge(:,new_edge(j)),local_Edge(1,3),2);
                    local_edge(1)   = old_edge(j);
                    local_Edge(:,1) = my_flip(Edge(:,old_edge(j)),local_Edge(1,2),2);
                end
            end
            n = n + 2;
        end
        % ----------------------- add one more internal point
        new_node = new_node + 1;
        xcoor=zeros(2,1); for j=1:2:7, xcoor = xcoor + Node(:,local_Edge(1,j)); end
        Node(:,nnode+new_node) = xcoor./4;

        % ----------------------- add four more internal edges
        Edge(1,CurEdge)   = local_Edge(2,1);
        Edge(2,CurEdge)   = nnode+new_node;
        Edge(1,CurEdge+1) = local_Edge(2,3);
        Edge(2,CurEdge+1) = nnode+new_node;
        Edge(1,CurEdge+2) = local_Edge(2,5);
        Edge(2,CurEdge+2) = nnode+new_node;
        Edge(1,CurEdge+3) = local_Edge(2,7);
        Edge(2,CurEdge+3) = nnode+new_node;

        % ----------- put the first quadrilateral at the place of the old one
        Face(:,iface) = [local_edge(1);CurEdge;CurEdge+3;local_edge(8)];
        % ----------------------- update faces
        Face(:,new_face)   = [local_edge(2);local_edge(3);CurEdge+1;CurEdge];
        Face(:,new_face+1) = [local_edge(4);local_edge(5);CurEdge+2;CurEdge+1];
        Face(:,new_face+2) = [local_edge(6);local_edge(7);CurEdge+3;CurEdge+2];

        %
        list  = [local_Edge(1,1);local_Edge(1,2);nnode+new_node; local_Edge(1,8)];
        [olist,fl] = order_face_nodes_quadTH(list,Node,1);
        Face_Node(1:4,iface,lvl+1)      = olist;

        list  = [local_Edge(1,2);local_Edge(1,3);local_Edge(1,4);nnode+new_node];
        [olist,fl] = order_face_nodes_quadTH(list,Node,1);
        Face_Node(1:4,new_face,  lvl+1) = olist;

        list  = [local_Edge(1,4);local_Edge(1,5);local_Edge(1,6);nnode+new_node];
        [olist,fl] = order_face_nodes_quadTH(list,Node,1);
        Face_Node(1:4,new_face+1,lvl+1) = olist;


        list = [local_Edge(1,6);local_Edge(1,7);local_Edge(1,8);nnode+new_node];
        [olist,fl] = order_face_nodes_quadTH(list,Node,1);
        Face_Node(1:4,new_face+2,lvl+1) = olist;

        Face_flag(new_face,  1) = Face_flag(iface,1);
        Face_flag(new_face+1,1) = Face_flag(iface,1);
        Face_flag(new_face+2,1) = Face_flag(iface,1);

        Face_thick(new_face,  1) = Face_thick(iface,1);
        Face_thick(new_face+1,1) = Face_thick(iface,1);
        Face_thick(new_face+2,1) = Face_thick(iface,1);

        Face_Parent(1:4,iface,lvl) = [iface; new_face; new_face+1; new_face+2];

        list  = local_Edge(1,[1,3,5,7]);
        [olist,fl] = order_face_nodes_quadTH(list,Node,1);
        Face_eorderP(1:4,iface)= olist;
        %    local_node_list4=order_face_nodes_quad(Node,olist,1);
        %    if local_node_list4~=olist, waitt('not ordered'),end

        mlist = [local_Edge(1,:) nnode+new_node];
        if fl==0,
            mlist=mlist([1,3,5,7,2,4,6,8,9]);
        else
            mlist=mlist([1,7,5,3,8,6,4,2,9]); %Main
        end

        %% Order the nodepoints as in the reference element, starting
        %  with the smallest coordinates in the left-lower corner
        [mxv,mxi] = find(Node(1,mlist)==min(Node(1,mlist)));
        [mvy,myi] = find(Node(2,mlist(mxi))==min(Node(2,mlist(mxi))));
        % The number of the point in 'mlist' is
        nm = mxi(myi);
        mlist_ref = [circshift(mlist(1:4),-(nm-1),2) circshift(mlist(5:8),-(nm-1),2) mlist(9)];

        %%
        Face_eorderD(1:9,iface)=mlist_ref;

        [l,~] = ashkan_reorder(Face_eorderD(:,iface), Node,9);
        Face_eorderD(:,iface) = l;
        [l,~] = ashkan_reorder(Face_eorderP(:,iface), Node,4);
        Face_eorderP(:,iface) = l;

        CurEdge  = CurEdge  + 4;
        new_face = new_face + 3;

    end   % ------------------------------ end loop over (coarse) faces
    nface_lvl(lvl+1)=size(Face,2);
    nnode_lvl(lvl+1)=size(Node,2);

    return
