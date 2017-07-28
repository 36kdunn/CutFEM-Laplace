function T = colorSquares(T,N,allCoords,segData)
%
% Labels squares 1 for exterior, 0 for intersected and -1 for interior in
% connectivity table T. Begins at each element marked 0 in the connectivity
% table and searches in the exterior direction to see if the neighboring
% square is incorrectly marked as interior (-1) when it should be
% exterior (1). The function searchNeighbors is called recursively until
% the current square is surrounded by 1s or 0s or if recurison limit is hit.
% Inputs:
%   T = connectivity table with column initialized to -1
%   N = number of intervals on one side of the square domain Omega.
%   allCoords = all coordinate points in G and where polygonal
%       approximation to G intersects the gridlines.
%   segData = starting points and number points in each curve segment.
% Outputs:
%   T = updated connectivity table

recLim = 7; % recursion limit to avoid filling stack
numSegs = length(segData); % number of curve segments

T(segData(:,3),10) = 0;

for sqInd = 1:N^2
    
    if T(sqInd,10) == 0
        
        % Find number of segments in current square
        segs = cell(sum(segData(:,3) == sqInd),1);
        ind = 1;
        
        % Collect all curve segments in one cell array
        for j = 1:numSegs
            if segData(j,3) == sqInd
                segs{ind} = allCoords(segData(j,1):segData(j,1)+segData(j,2)-1,:);
                ind = ind+1;
            end
        end
        
        % Store current segment as S and define necessary quantities
        S = segs{1};
        len = length(S); % number of points in curve segment
        % Midpoint of linear segment connecting endpoints
        mdpt = [S(len,1)+S(1,1),+(S(len,2)+S(1,2))]/2;
        % Outward normal to linear segment connecting endpoints
        n = [S(len,2)-S(1,2),-(S(len,1)-S(1,1))];
        
        % Define x and y search directions based on signs of outward normal
        xDir = sign(n(1,1));
        yDir = sign(n(1,2));
        
        for j = 1:length(segs)
            
            Stemp = segs{j};
            len = length(Stemp);
            mdpttemp = [S(len,1)+S(1,1),+(S(len,2)+S(1,2))]/2;
            
            if xDir*mdpttemp(1,1) > xDir*mdpt(1,1)
                xDir = 0;
            end
            
            if yDir*mdpttemp(1,2) > yDir*mdpt(1,2)
                yDir = 0;
            end
            
        end
        
        % Check neighbors in appropriate directions
        if yDir == 1 && sqInd + N <= N^2 && T(sqInd+N,10)~=0
            T = searchNeighbors(T,N,sqInd+N,recLim);
        end
        if yDir == -1 && sqInd - N > 0 && T(sqInd-N,10)~=0
            T = searchNeighbors(T,N,sqInd-N,recLim);
        end
        if xDir == 1 && sqInd + 1 <= N^2 && T(sqInd+1,10)~=0
            T = searchNeighbors(T,N,sqInd+1,recLim);
        end
        if xDir == -1 && sqInd - 1 > 0 && T(sqInd-1,10)~=0
            T = searchNeighbors(T,N,sqInd-1,recLim);
        end
        
    end % switch
    
end

% Search squares starting with the edges of Omega.  This is to make sure we
% don't miss anything when we simply emanate from the interface.
for rowInd = 1:N
    
    flagLeft = 1;
    flagRight = 1;
    flagUp = 1;
    flagDown = 1;
    
    for colInd = 1:ceil(N/2)
        
        leftInd = colInd + N*(rowInd-1);
        if T(leftInd,10) == -1 && flagLeft == 1
            T(leftInd,10) = 1;
        else
            if T(leftInd,10) == 0
                flagLeft = 0;
            end
        end
        
        rightInd = (N-(colInd - 1)) + N*(rowInd-1);
        if T(rightInd,10) == -1 && flagRight == 1
            T(rightInd,10) = 1;
        else
            if T(rightInd,10) == 0
                flagRight = 0;
            end
        end
        
        upInd = rowInd + N*(colInd-1);
        if T(upInd,10) == -1 && flagUp == 1
            T(upInd,10) = 1;
        else
            if T(upInd,10) == 0
                flagUp = 0;
            end
        end
        
        downInd = N^2 - (rowInd-1) - N*(colInd-1);
        if T(downInd,10) == -1 && flagDown == 1
            T(downInd,10) = 1;
        else
            if T(downInd,10) == 0
                flagDown = 0;
            end
        end
        
    end
    
end

end

function T = searchNeighbors(T,N,index,recLim)
%
% Recursive function checking neighbors if appropriate unti recursion limit

if recLim == 0
    return
end

if T(index,10) == 0
    return
end

if T(index,10) == -1
    T(index,10) = 1;
end

% search down
if (index - N) > 0
    % If neighbor is marked as interior, it must be by default
    if T(index - N,10) == -1
        T = searchNeighbors(T,N,index - N,recLim-1);
    end
end

% search up
if (index + N) <= N^2
    % If neighbor is marked as interior, it must be by default
    if T(index + N,10) == -1
        T = searchNeighbors(T,N,index + N,recLim-1);
    end
end

% search left
if (index - 1) >0
    % If neighbor is marked as interior, it must be by default
    if T(index - 1,10) == -1
        T = searchNeighbors(T,N,index - 1,recLim-1);
    end
end

% search right
if (index + 1) < N^2+1
    % If neighbor is marked as interior, it must be by default
    if T(index + 1,10) == -1
        T = searchNeighbors(T,N,index + 1,recLim-1);
    end
end

end

