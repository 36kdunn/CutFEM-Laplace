function [polyCoords,polyIndices] = findPolygons(V,N,segs,segInd,polyIndices,polyCoords,Tol)
%
% This function searches over the points provided in V, segs to label the
% points in a counter-clockwise direction, separating them into the
% polygons that are make up the area of the square formed by V.
%
% Inputs:
%   V: (four) vertices of square, labelled in a counter-clockwise manner???
%   segs: vertices of the segment, labelled counter-clockwise
%   segInd: index of the last found segment
%   polyIndices: [indCW, indACW, sign] where indCW is the index of
%             clockwise-oriented polygons, indACW is the index of the
%             counterclockwise-oriented polygons, and sign is 1 if
%             counterclockwise and -1 if clockwise
%   polyCoords: cell array of indices for polygons; 1-4 are clockwise
%             oriented and 5-8 are counterclockwise-oriented
%   Tol: error tolerance for comparison, i.e., x=y if abs(x-y)<Tol
%
% Outputs:
%   polyCoords: output for recursive calls to store additional points or
%             additional polygons
%   polyIndices: output for recursive calls to know which is the next
%             clockwise or counterclockwise index and the orientation of
%             the next call

% Store index values for this current iteration. Important: these will not
% be changed by recursive calls.
sgn = polyIndices(3); % -1 if this polynomial is oriented clockwise, 1 if counterclockwise
polyInd = polyIndices(1+(1+sgn)/2); % 1 if clockwise, 2 if counterclockwise

% Insert curve segment as first part of polygon
polyCoords{polyInd} = [polyCoords{polyInd};segs{segInd}];

if polyInd == 5
    % This is the first iteration of the loop and must be called on the
    % same segment in both counterclockwise and clockwise directions
    PI_temp = [1,5,-1]; % Same as initial condition, but clockwise
    [polyCoords,polyIndices] = findPolygons(V,N,segs,1,PI_temp,polyCoords,Tol);
    polyIndices(3) = -polyIndices(3); % reverse direction to avoid messing up recursive calls
end

% Save x and y values where loop should break if encountered
endx = segs{segInd}(1,1);
endy = segs{segInd}(1,2);

% Initialize condition and begin loop
break_flag = 0;
iter = 1;
iterMax = 100;
while break_flag == 0 && iter<iterMax
    
    % Initialize skip condition for each iteration. If flag==1, then a
    % point/segment has been added to the polygon and we start back at
    % the beginning of the loop.
    flag = 0;
    
    % Save endpoint of segment
    x = polyCoords{polyInd}(end,1);
    y = polyCoords{polyInd}(end,2);
    
    % Outward normal of counterclockwise oriented segment is [dy,-dx]
    n = [y-polyCoords{polyInd}(end-1,2), polyCoords{polyInd}(end-1,1)-x];
    % If any component is close enough to zero, change it to zero to avoid
    % sign errors in the normal direction
    n(abs(n)<Tol) = 0; 
    
    % Search counterclockwise from this point
    if abs(y - round(y*N)/N) < Tol
        % Current endpoint is on top or bottom of square
        % Begin by checking for curve segments along the top/bottom of the square
        [polyCoords,polyIndices,flag,break_flag] = ...
            checkSegsTopBot(V,N,polyCoords,segs,polyIndices,x,y,n,polyInd,sgn,endx,endy,Tol);
        
        % If we have not found a segment, we check the top/bottom corners of the square
        if flag == 0 && break_flag == 0
            [polyCoords,flag] = checkCornersTopBot(V,polyCoords,polyInd,x,y,n,sgn,Tol);
        end
    end
    
    % Search counterclockwise from this point
    if abs(x - round(x*N)/N) < Tol && flag == 0 && break_flag == 0
        % Current endpoint is on left or right of square
        
        % First check for curve segments along the left/right sides of the square
        [polyCoords,polyIndices,flag,break_flag] = checkSegsLR(V,N,polyCoords,segs,polyIndices,x,y,n,polyInd,sgn,endx,endy,Tol);
        
        % If we have not found a segment, we check the left/right corners of the square
        if flag == 0 && break_flag == 0
            polyCoords = checkCornersLR(V,polyCoords,polyInd,x,y,n,sgn,Tol);
        end
    end
    iter = iter + 1;
end

if iter>=iterMax
    error(['Find polygon loop hung up in square with bottom left corner x=',num2str(V(1,1)),'   y=',num2str(V(1,2))]);
end

end


function [polyCoords,polyIndices,flag,break_flag] = checkSegsTopBot(V,N,polyCoords,segs,polyIndices,x,y,n,polyInd,sgn,endx,endy,Tol)
%
% Note on inputs:
%   polyIndices is inherited to pass to recursive calls, but polyInd and
%   sgn are also inherited since these are the values of the current
%   polygon and any future additions must be appended as such.

% Initialize values and begin loop
break_flag = 0; % turns to 1 if the first point of segment is found
flag = 0; % turns to 1 if we find/add a segment
idx = 0; % stores the index of the first new segment to be added
currentX = x + 1; % initialized as a value outside the domain
for i = 1:length(segs)
    % Looks to see if
    %   (1) current endpoint is in direction opposite x comp of normal
    %   (2) the segment is on the same edge of the square
    if abs(segs{i}(1,2) - y) < Tol && sign(x-segs{i}(1,1)) == sgn*sign(n(1,1))
        
        % If this is the closest endpoint observed so far
        if abs(segs{i}(1,1)-x) < abs(currentX - x)
            % Update idx to current index and currentX to x value
            idx = i;
            currentX = segs{i}(1,1);
        end
    end
end

if idx > 0
    % If endpoint was found
    if abs(segs{idx}(1,1)-endx) + abs(segs{idx}(1,2)-endy) < 2*Tol
        % Check to see if this is the termination endpoint
        break_flag = 1;
    else
        % Otherwise, add segment to coordinates of polygon
        flag = 1;
        polyCoords{polyInd} = [polyCoords{polyInd};segs{idx}];
        
        % Reverse direction
        polyIndices(3) = -polyIndices(3);
        % Update the index in the reverse direction
        polyIndices(1+(1+polyIndices(3))/2) = polyIndices(1+(1+polyIndices(3))/2) + 1;
        % Make recursive call from endpoint of segment in opposite direction
        [polyCoords,polyIndices] = findPolygons(V,N,segs,idx,polyIndices,polyCoords,Tol);
        % Reverse direction again to avoid messing up previous recursive calls
        polyIndices(3) = -polyIndices(3);
    end
end

end


function [polyCoords,polyIndices,flag,break_flag] = checkSegsLR(V,N,polyCoords,segs,polyIndices,x,y,n,polyInd,sgn,endx,endy,Tol)
%
% Note on inputs:
%   polyIndices is inherited to pass to recursive calls, but polyInd and
%   sgn are also inherited since these are the values of the current
%   polygon and any future additions must be appended as such.

% Initialize values and begin loop
break_flag = 0; % turns to 1 if the first point of segment is found
flag = 0; % turns to 1 if we find/add a segment
idx = 0; % stores the index of the first new segment to be added
currentY = y + 1; % initialized as a value outside the domain
for i = 1:length(segs)
    % Looks to see if
    %   (1) current endpoint is in direction opposite y comp of normal
    %   (2) the segment is on the same edge of the square
    if abs(segs{i}(1,1) - x) < Tol && sign(y-segs{i}(1,2)) == sgn*sign(n(1,2))
        
        % If this is the closest endpoint observed so far
        if abs(segs{i}(1,2)-y) < abs(currentY - y)
            % Update idx to current index and currentY to y value
            idx = i;
            currentY = segs{i}(1,2);
        end
    end
    
end

% If endpoint was found
if idx > 0
    % Check to see if this is the termination endpoint
    if abs(segs{idx}(1,1)-endx) + abs(segs{idx}(1,2)-endy) < 2*Tol
        break_flag = 1;
    else
        % Otherwise, add segment to coordinates of polygon
        flag = 1;
        polyCoords{polyInd} = [polyCoords{polyInd};segs{idx}];
        
        % Reverse direction
        polyIndices(3) = -polyIndices(3); % reverse direction
        % Update the index in the reverse direction
        polyIndices(1+(1+polyIndices(3))/2) = polyIndices(1+(1+polyIndices(3))/2) + 1;
        % Make recursive call from endpoint of segment in opposite direction
        [polyCoords,polyIndices] = findPolygons(V,N,segs,idx,polyIndices,polyCoords,Tol);
        % Reverse direction again to avoid messing up previous recursive calls
        polyIndices(3) = -polyIndices(3);
    end
end

end


function [polyCoords,flag] = checkCornersTopBot(V,polyCoords,polyInd,x,y,n,sgn,Tol)

% Initialize to zero
flag = 0;

% Loop over vertices
for i = 1:numel(V)/2
    % Looks to see if
    %   (1) current vertex is in direction opposite x comp of normal
    %   (2) the segment is on the same edge of the square
    if abs(V(i,2)- y) < Tol && sign(x-V(i,1)) == sgn*sign(n(1,1))
        % If current vertex is not the same as the endpoint of a segment
        if V(i,1) ~= polyCoords{polyInd}(end,1) || V(i,2) ~= polyCoords{polyInd}(end,2)
            % Add corner to polygon, change flag, and break loop
            if abs(V(i,1)- x) > Tol % Make sure corner is not equal to a previous curve seg endpoint
                polyCoords{polyInd} = [polyCoords{polyInd};V(i,:)];
                flag = 1;
                break
            end
        end
    end
end




end


function polyCoords = checkCornersLR(V,polyCoords,polyInd,x,y,n,sgn,Tol)

% Loop over vertices
for i = 1:numel(V)/2
    % Looks to see if
    %   (1) current vertex is in direction opposite x comp of normal
    %   (2) the segment is on the same edge of the square
    if abs(V(i,1)- x) < Tol && sign(y-V(i,2)) == sgn*sign(n(1,2))
        % If current vertex is not the same as the endpoint of a segment
        if V(i,1) ~= polyCoords{polyInd}(end,1) || V(i,2) ~= polyCoords{polyInd}(end,2)
            % Add corner to polygon, change flag, and break loop
            if abs(V(i,2)- y) > Tol % Make sure corner is not equal to a previous curve seg endpoint
                polyCoords{polyInd} = [polyCoords{polyInd};V(i,:)];
                break
            end
        end
    end
end

end