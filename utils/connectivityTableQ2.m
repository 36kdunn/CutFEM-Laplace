function [T,V] = connectivityTableQ2(N)
%
% Construct connectivity table for Q2 elements.  The numbering of the nodes
% is as follows:
%   3 - 6 - 9
%   |   |   |
%   2 - 5 - 8
%   |   |   |
%   1 - 4 - 7
% Global numbering of the nodes begins at 1 in the bottom left corner of
% the domain and progresses right until the end is reached, then up to the
% left-most point on the above line.
%
% Inputs:
%   N = number of intervals on one side of the square domain Omega.
% Outputs:
%   T = table of global vertex numbers.  Each row corresponds to an
%      element where the columns are numbered as in the diagram above.
%      The 10th column is -1 for interior elements, 1 for exterior elements,
%      and 0 for elements intersected by Gamma; initialized as -1 for all elements.
%   V = table of coordinates of the global nodes.  The ith row contains
%      the (x,y) coordinates of node i in the global numbering.

% Global numbering of bottom-left corner of every square
x = zeros(N^2,1);
x(1:N,1) = 1:2:(2*N-1);

% Number nodes for corners of rectangles
for i = 0 : N-1
    x(i*N + 1 : (i+1)*N) = (1:2:(2*N-1)) + i*2*(2*N+1);
end

% Numbering of nodes in each element
T = zeros(N^2,10);
T(:,1:3) = [x, x+(2*N+1), x+2*(2*N+1)]; % bot left, middle left, top left
T(:,4:6) = T(:,1:3) + 1; % bot middle, middle middle, top middle
T(:,7:9) = T(:,1:3) + 2; % bot right, middle right, top right
T(:,10) = -1;

% Coordinates of nodes
V = [kron(ones(1,2*N+1),0:1/(2*N):1)',kron(0:1/(2*N):1,ones(1,2*N+1))'];


end