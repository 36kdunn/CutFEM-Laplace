function [allCoords,segData] = findInterfaceIntersections(G,N)
%
% Integrates bulk basis functions over elements intersected by the interface.
% This function is defined for interior and exterior elements.
% Line segments are integrated counterclockwise and coordinates must be given
% counerclockwise.
%
% Inputs:
%    G = interface by coordinates
%    N = number of divisions on x- and y-axis
%
% Outputs:
%    allCoords = coordinates where Gamma intersects grid lines drawn by
%       elements and original points in G, ordered counterclockwise.
%    segData = data necessary to construct curve segments of Gamma between
%       grid line intersections. Column 1 contains starting point of each
%       curve segment corresponding to the rows of allCoords and column 2
%       contains the number of points on this segment.

% Initialize loop parameters
ind = 1;
coordsTemp = G(1,:);
segData = [];

% Find initial point of curve segment, i.e., coordinates of point where
% curve intersects grid line in first element.
% Also output index of last point in G inside element.
[coordsTemp,Gmax] = findSegCoords(G,N,1,coordsTemp,-1);
if norm(coordsTemp(1,:) - coordsTemp(2,:)) < eps
    coordsTemp = coordsTemp(2:end,:);
end

allCoords = [];
curveCounter = 1;

% Modify last curve index to pass to while loop and findCoordinates function
if Gmax == 1
    Gmax = length(G);
else
    Gmax = Gmax - 1;
end

iterMax = N^2;
iter = 1;

% Loop over elements.  If ind == Gmax, then the loop has fully traced G
while ind < Gmax
    
    [coordsTemp,ind] = findSegCoords(G,N,ind,coordsTemp,1);
    if norm(coordsTemp(end,:) - coordsTemp(end-1,:)) < eps
        coordsTemp = coordsTemp(1:end-1,:);
    end
    if numel(coordsTemp) > 4 && norm(coordsTemp(1,:) - coordsTemp(2,:)) < eps
        coordsTemp = coordsTemp(2:end,:);
    end
    
    if ~isempty(allCoords)
        if norm(coordsTemp(1,:) - allCoords(end,:)) < eps
            [numSegs,~] = size(coordsTemp); % order matters here
            midpt = mean(coordsTemp,1);
            coordsTemp = coordsTemp(2:end,:);
        else
            [numSegs,~] = size(coordsTemp);
            midpt = mean(coordsTemp,1);
        end
        [curveInd,~] = size(allCoords);
    else
        curveInd = 1;
        [numSegs,~] = size(coordsTemp);
        midpt = mean(coordsTemp,1);
    end
    
    
    % Compute global index of x and y in grid
    xind = floor(midpt(1,1)*N + eps);
    yind = floor(midpt(1,2)*N + eps);
    
    segGlobalInd = (xind+1) + N*yind;
    allCoords = [allCoords;coordsTemp];
    segData = [segData; [curveInd, numSegs, segGlobalInd]];
    
    % Prepare for next iteration (if there is one)
    coordsTemp = allCoords(end,:);
    curveCounter = curveCounter+1;
    iter = iter+1;
    
    if iter > iterMax
        error(['findInterfaceIntersections hung up at x=',num2str(G(ind,1)),', y=',num2str(G(ind,2))]);
    end
end

% Delete errors
flag = 0;
ind = 1;
len = length(segData);
while ind < len
    if segData(ind,2) == 1
        if ind < len
            if segData(ind,1) == segData(ind+1,1)
                flag = 1;
            end
        end
        if ind > 1
            if segData(ind,1) == segData(ind-1,1)
                flag = 1;
            end
        end
    end
    if flag == 1
        segData = [segData(1:ind-1,:);segData(ind+1:len,:)];
    end
    len = length(segData);
    ind = ind+1;
    flag = 0;
end

end
