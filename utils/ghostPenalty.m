function A = ghostPenalty(T,N,A,firstInd,localIndsU,gam2,gam3,M,intext)

% gam3=0;

% AP = A;

for i = 1:N^2
    
    condVert = 0;
    condHoriz = 0;
    
    if mod(i,N) > 0 && i+1 <= N^2
        right = i + 1;
        if intext == 1
            condVert = (T(i,10) == 0 && T(right,10) >= 0) || (T(i,10) == 1 && T(right,10) == 0);
%             condVert = (T(i,10) >= 0 && T(right,10) >= 0);
        else
            condVert = (T(i,10) == 0 && T(right,10) <= 0) || (T(i,10) == -1 && T(right,10) == 0);
%             condVert = (T(i,10) <= 0 && T(right,10) <= 0);
        end
    end
    
    if i + N <= N^2
        up = i + N;
        if intext == 1
            condHoriz = (T(i,10) == 0 && T(up,10) >= 0) || (T(i,10) == 1 && T(up,10) == 0);
%             condHoriz = (T(i,10) >= 0 && T(up,10) >= 0);
        else
            condHoriz = (T(i,10) == 0 && T(up,10) <= 0) || (T(i,10) == -1 && T(up,10) == 0);
%             condHoriz = (T(i,10) <= 0 && T(up,10) <= 0);
        end
        
    end
    
    % % % VERTICAL EDGES
    if condVert
        colGlobalInds = [T(i,1:9),T(right,4:9)];
        
        inds = firstInd + localIndsU(colGlobalInds);
        A(inds,inds) = A(inds,inds) + gam2*M.ghostvert;
        A(inds,inds) = A(inds,inds) + gam3*M.ghostvert2ndOrder;
    end
    
    % % % HORIZONTAL EDGES
    if  condHoriz
        colGlobalInds = [T(i,1:3), T(up,2:3), T(i,4:6), T(up,5:6), T(i,7:9), T(up,8:9)];
        
        inds = firstInd + localIndsU(colGlobalInds);
        A(inds,inds) = A(inds,inds) + gam2*M.ghosthoriz;
        A(inds,inds) = A(inds,inds) + gam3*M.ghosthoriz2ndOrder;
    end
    
end

end

