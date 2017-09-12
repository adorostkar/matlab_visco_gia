function [X,Y, id] = reorder_nodes(Node)
% Create the matrices of the correct size to fill.
% X is the y coordinates of the Nodes
% Y is the y coordinates of the Nodes
% id is the index of the corresponding element in Node list. This will be
% used to reorder the solution for the mesh plot
%
%   I M P O R T A N T    N O T I C E
% This routin relies on the nodes being in the same place (on a line) in each y
%  line. This should be used exactly after the mesh is constructed.
% if the mesh moves, this routin will not work anymore. So
% it has to be used once first to have the locations and then the saved locations be used in the rest of the application

% X = zeros(2*Ny-1,2*Nx-1);
% Y = zeros(2*Ny-1,2*Nx-1);
% id= zeros(2*Ny-1,2*Nx-1);
%
% for i = 1:Ny
% id(2*i-1, 1:2:end) = (i-1)*Nx+1:i*Nx;
% end
%
% for i = 1:Ny
% id(2*i-1, 2:2:end) = Nx*Ny + (i-1)*(2*Nx-1) + 1: Nx*Ny + (i-1)*(2*Nx-1) + Nx - 1;
% end
%
% for i = 1:Ny-1
% id(2*i, 1:2:end) = Nx*Ny + i*(Nx-1) + (i-1)*Nx + 1: Nx*Ny + i*(Nx-1) + (i-1)*Nx + Nx;
% end
%
% for i = 1:Ny-1
% id(2*i, 2:2:end) = (2*Nx-1)*Ny + Nx*(Ny-1) +(i-1)*(Nx-1)+ 1: (2*Nx-1)*Ny + Nx*(Ny-1) + i*(Nx-1);
% end
%
% X(:) = Node(1,id(:));
% Y(:) = Node(2,id(:));
tol = 1e-16;
xx = [];
yy = [];
Ixx = [];
Iyy = [];
while (length(Node(1,:)) ~= length(xx)*length(yy) && tol < 1)
    [xx, ~,Ixx] = uniquetol(Node(1,:), 1e-15);
    [yy, ~,Iyy] = uniquetol(Node(2,:), 1e-15);
    tol = tol*10;
end

% Sort y direction so that the largest value is last which produces indecis
% in the natural vector direction.
%id(1,:) is where y is the largest and all x values
[~, I] = sort(yy,'descend');
yy = yy(I);
Iyy = I(Iyy);

if tol == 1
    error('Unable to find the correct number of unique coordinates');
end

[X, Y] = meshgrid(xx, yy);
id = zeros(size(X));



for i = 1:size(Node,2)
    ii = Iyy(i);
    jj = Ixx(i);
    id(ii,jj) = i;
    % Z(ii,jj) = vv;
end

end
