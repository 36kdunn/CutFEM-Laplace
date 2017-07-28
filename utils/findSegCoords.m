function [coords,Gind] = findSegCoords(G,N,Gind,coords,flag)
%
% Find the points of intersections of the interface with the gridlines by
% recursively calling searchInterface. This will search counterclockwise
% along interface if flag==1 and clockwise if flag==-1.
% Inputs:
%   G = original points along interface Gamma
%   N = number of intervals on one side of the square domain Omega.
%   Gind = index of current point in G
%   coords = current list of (x,y) coordinates along Gamma, including
%       intersections with grid.
%   flag = 1 to search counterclockwise, -1 to search clockwise.
% Outputs:
%   coords = list of coordinates along Gamma including intersections with
%       grid, now updated with latest segment.
%   Gind = index of next starting in G.

% Define step size of grid (in x and y)
h = 1/N;

updatedInd = Gind + flag;
if updatedInd <= 0
    updatedInd = length(G) - 1;
end


% Define temporary endpoints of line segment
% Starting point
x1 = coords(end,1);
y1 = coords(end,2);
% In direction of
x2 = G(updatedInd,1);
y2 = G(updatedInd,2);

% Find location in grid
% NOTE: TOLERANCE
xloc = floor(x1/h + N*eps);
yloc = floor(y1/h + N*eps);


% Find x bounds of current rectangle. Determines which direction to go in
% if (x1,y1) is on an x grid line
if x1 == xloc*h && flag*sign(x2-x1) < 0
    Mx = xloc*h; % max
    mx = (xloc-1)*h; % min
else
    Mx = (xloc+1)*h; % max
    mx = xloc*h; % min
end

% Find y bounds of current rectangle. Determines which direction to go in
% if (x1,y1) is on a y grid line
if y1 == yloc*h && flag*sign(y2-y1) < 0
    My = yloc*h; % max
    my = (yloc-1)*h; % min
else
    My = (yloc+1)*h; % max
    my = yloc*h; % min
end

% Package intervals and pass to recursive search function
intervals = [[mx,Mx];[my,My]];
% while min(abs(coords(1,1)-intervals(1,:))) > 0 && min(abs(coords(1,2)-intervals(2,:))) > 0
[coords,Gind] = searchInterface(G,Gind,coords,intervals,flag);
%     Gind = Gind - 1;
% end
% Gind = Gind + 1;

end

function [coords,Gind] = searchInterface(G,Gind,coords,intervals,flag)

% If counterclockwise, then most recent point on curve will be stored at
% end of the coordinate vector.  If clockwise, at the beginning.
if flag == 1
    xc = coords(end,1);
    yc = coords(end,2);
else
    xc = coords(1,1);
    yc = coords(1,2);
end

% Update index in clockwise/counterclockwise direction, if flag is -1/1
updatedInd = Gind + flag;
if updatedInd <= 0
    updatedInd = length(G) - 1;
end
if updatedInd > length(G)
    updatedInd = 2; % avoid index out of bounds error on final iteration
end

% Store new endpoint (direction) of line segment
xn = G(updatedInd,1);
yn = G(updatedInd,2);

% Name interval bounds
mx = intervals(1,1); % min in x direction
Mx = intervals(1,2); % max in x direction
my = intervals(2,1); % min in y direction
My = intervals(2,2); % max in y direction

[len,~]=size(coords);

% If the saved point is located in the element/on the boundary of the
% element, concatenate coordinate matrix appropriately and recursively call
% searchInterfacce function.
% Else, find intersection of line segment from (xc,yc) in the direction of
% (xn,yn) and boundary of element
if (xn<= Mx && xn >= mx) && (yn <= My && yn >= my)
    if flag == 1
        coords = [coords;[xn,yn]];
    else
        coords = [[xn,yn];coords];
    end
    [coords,Gind] = searchInterface(G,updatedInd,coords,intervals,flag);
    
else
    % If this is not the original point of the segment and the point lies
    % on a grid line, but the grid line is not included in the current
    % intervals, then don't calculate/store new point and return.
    if len > 1 && ((xc == mx || xc == Mx) || (yc == my || yc == My))
        return
    end
    
    % Compute angle between positive x axis and line segment
    theta = atan2((yn-yc),(xn-xc)); % compute angle made by line segment
    
    % Assume segment intersects left or right side of cell.  Compute
    % (xp,yp) and ds (length of segment) accordingly.
    if sign(xn-xc) < 0
        xp = mx;
    else
        xp = Mx;
    end
    ds = abs((xc-xp)/cos(theta)); % calculate distance between (xp,yp) and (xn,yn)
    
    if xn == xc
        % If x1 == x2 then cos(theta) = 0 and bad stuff happens with ds and
        % we can't compute yp.  Instead, yp will be defined to be boundary
        % of interval in direction of segment.
        if sign(yn-yc) < 0
            yp = my;
        else
            yp = My;
        end
    else
        yp = yc + ds*sin(theta); % define "next" y point in sequence
    end
    
    % Assume segment intersects top or bottom of cell. Compute y coordinate
    % yp1 and ds1 (length of prospective segment) accordingly.  The x
    % coordinate will be defined later if necessary.
    if sign(yn-yc) < 0
        yp1 = my;
    else
        yp1 = My;
    end
    ds1 = abs((yc-yp1)/sin(theta));
    
    % Both line segments are in the same direction.  Shortest line segment
    % must intersect respective boundary first and is chosen.  Segment with
    % length ds is chosen by default.  If necessary, switch to segment with
    % length ds1.
    if ds1 < ds || ds == 0
        % if segment intersects top or bottom of cell, redefine variables
        if yc == yn
            % if y1 == y2 then sin(theta) = 0 and bad stuff happens with ds1
            if sign(xn-xc) < 0
                xp = mx;
            else
                xp = Mx;
            end
        else
            xp = xc + ds1*cos(theta);
        end
        yp = yp1;
    end
    
    % Place new coordinates in correct place
    if flag == 1
        coords = [coords;[xp,yp]];
    else
        coords = [[xp,yp];coords];
    end
end

end