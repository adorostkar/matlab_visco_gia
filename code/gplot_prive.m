function [Xout,Yout]=gplot_prive(A,xy,color,lc)
%GPLOT Plot graph, as in "graph theory".
%   GPLOT(A,xy) plots the graph specified by A and xy. A graph, G, is
%   a set of nodes numbered from 1 to n, and a set of connections, or
%   edges, between them.  
%
%   In order to plot G, two matrices are needed. The adjacency matrix,
%   A, has a(i,j) nonzero if and only if node i is connected to node
%   j.  The coordinates array, xy, is an n-by-2 matrix with the
%   position for node i in the i-th row, xy(i,:) = [x(i) y(i)].
%   
%   GPLOT(A,xy,LineSpec) uses line type and color specified in the
%   string LineSpec. See PLOT for possibilities.
%
%   [X,Y] = GPLOT(A,xy) returns the NaN-punctuated vectors
%   X and Y without actually generating a plot. These vectors
%   can be used to generate the plot at a later time if desired.
%   
%   See also SPY, TREEPLOT.

%   John Gilbert, 1991.
%   Modified 1-21-91, LS; 2-28-92, 6-16-92 CBM.
%   Copyright (c) 1984-96 by The MathWorks, Inc.
%   $Revision: 5.6 $  $Date: 1996/10/24 19:01:23 $

[i,j] = find(A);
[ignore, p] = sort(max(i,j));
i = i(p);
j = j(p);

% Create a long, NaN-separated list of line segments,
% rather than individual segments.

X = [ xy(i,1) xy(j,1) repmat(NaN,size(i))]';
Y = [ xy(i,2) xy(j,2) repmat(NaN,size(i))]';
X = X(:);
Y = Y(:);

if nargout==0,
    if nargin<3,
        H=plot(X, Y);
        set(H,'LineWidth',3,'Color','b')
    elseif nargin<4,
        H=plot(X, Y);
        set(H,'LineWidth',3,'Color',color)
    else
        H=plot(X, Y, lc);
        set(H,'LineWidth',3,'Color',color)
    end
else
    Xout = X;
    Yout = Y;
end
