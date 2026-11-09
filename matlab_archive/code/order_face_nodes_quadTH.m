% Order the edges/nodes in a given face so that the nodes
% form a closed path and are ordered
% clockwise         if clcws = 1, or
% counterclockwise  if clcws = 0
%
% ------------------- make an oriented local copy of the nodes
% ------------------- use a relation Face_Node(:,iface,lvl)


function [local_node_list,fl]=...
        order_face_nodes_quadTH(local_node_list,Node,clcws)
    
    fl = 0;
    %%% --- clock / counterclock orientation
    vec_prod = ( (Node(1,local_node_list(2))-Node(1,local_node_list(1)))*...
        (Node(2,local_node_list(3))-Node(2,local_node_list(1))) ) - ...
        ( (Node(1,local_node_list(3))-Node(1,local_node_list(1)))*...
        (Node(2,local_node_list(2))-Node(2,local_node_list(1))) );
    
    if vec_prod > 0, % current ordering is counterclock-wise
        if clcws == 1
            qq=local_node_list(2);local_node_list(2)=local_node_list(4);
            local_node_list(4)=qq;
            fl = 1;
        end
    end
    if vec_prod < 0, % current ordering is clock-wise
        if clcws == 0
            qq=local_node_list(2);local_node_list(2)=local_node_list(4);
            local_node_list(4)=qq;
        end
    end
    % additional shift to have the smallest node number
    % (which corresponds to the coarse node) first
    
    [iv,ii]=min(local_node_list);
    local_node_list = circshift(local_node_list,[-(ii-1),0]); % shift
    
    return
    %
