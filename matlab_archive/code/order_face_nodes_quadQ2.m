% Order the nodes in a given face so that the nodes
% form a closed path and are ordered
% Presummed ordering of the nodepoints in the reference f.e.:
% (-1,1)      (1,1)
%
%    2    6     3
%     ---------
%    |	  |    |
%    |	  |    |
%   5|____9____|7      FUN(9)
%    |	  |    |
%    |	  |    |
%   1|____|____|4
%         8
% The first 4 points in local_node_list are already ordered.
% The list is sorted in an increasing order,
% i.e., point nr. 9 is in its place.
% Points 5-8 are numbered via Node_Edge.
%


function [local_node_list9]=...
        order_face_nodes_quadQ2_new(Node,local_node_list4,local_node_list90)
    
    local_node_list9 = [local_node_list4;0;0;0;0;0];
    v = [1,2,3,4,1];
    Coor_n=Node(:,local_node_list90(5:9));
    for k=1:4,
        Coor_o=Node(:,[local_node_list4(v(k)),local_node_list4(v(k+1))]);
        Mid(:,k) = 0.5*sum(Coor_o');
        p=find(((abs(Coor_n(1,:)-Mid(1,k))<1e-6))&((abs(Coor_n(2,:)-Mid(2,k))<1e-6)));
        local_node_list9(4+k) = local_node_list90(p+4);
    end
    local_node_list9(9) = max(local_node_list90);
    
    return
    %
