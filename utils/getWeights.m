function [Wext,Wint] = getWeights(T,V,N,allCoords,segData,extLocalInds,intLocalInds,Tol)
%

Nue = sum(extLocalInds>0);
Wext = zeros(Nue,1);

Nui = sum(intLocalInds>0);
Wint = zeros(Nui,1);

integrals = [1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36]'*(1/N)^2;

numSegs = length(segData);

[numSq,~] = size(T);

for sqInd = 1:numSq
    
    switch T(sqInd,10)
        
        case 1 % Exterior square
            extInds = extLocalInds(T(sqInd,1:9),1);
            Wext(extInds,1) = Wext(extInds,1) + integrals;
            
        case -1 % Interior square
            intInds = intLocalInds(T(sqInd,1:9),1);
            Wint(intInds,1) = Wint(intInds,1)  + integrals;
            
        case 0 % Interface square
            
            segs = cell(sum(segData(:,3) == sqInd),1);
            ind = 1;
            
            % Collect all curve segments in one cell array to pass to findPolygons
            for j = 1:numSegs
                if segData(j,3) == sqInd
                    segs{ind} = allCoords(segData(j,1):segData(j,1)+segData(j,2)-1,:); % maybe coordsInds(j,2)-2 ???
                    ind = ind+1;
                end
            end
            
            polyCoords = cell(8,1);
            polyInd = [0, 5, 1];
            segInd = 1;
            [polyCoords] = findPolygons(V(T(sqInd,[1,3,7,9]),:),N,segs,segInd,polyInd,polyCoords,Tol);
            
            
            % Translate RHS functions to current square
            x0 = V(T(sqInd,1),1);
            y0 = V(T(sqInd,1),2);
            fWTrans = @(x,y) fW(x,y,x0,y0,N);
            
            [Iext,Iint] = integrateCase0(fWTrans,polyCoords,Tol);
            
            % Get local exterior indices
            extInds = extLocalInds(T(sqInd,1:9),1);
            Wext(extInds,1) = Wext(extInds,1) + Iext;
            
            % Get local interior indices
            intInds = intLocalInds(T(sqInd,1:9),1);
            Wint(intInds,1) = Wint(intInds,1) + Iint;

    end % switch
    
end

end

function I = fW(x,y,x0,y0,N)

I = basisFuncsVec(x,y,x0,y0,N);

end

function [Iext,Iint] = integrateCase0(f,polyCoords,Tol)

Iext = zeros(9,1);
Iint = Iext;


for k = 1:4
    if ~isempty(polyCoords{k})
        % First four cells in polyCoords are ordered clockwise
        % which means they belong to the exterior
        verts = zeros(size(polyCoords{k}));
        [N,~] = size(verts);
        for l = 1:N
            verts(l,:) = polyCoords{k}(N-l+1,:);
        end
        Iext = Iext + intOverPoly(f,verts,Tol);
    end
    
    if ~isempty(polyCoords{k+4})
        % Last four cells in polyCoords are ordered counterclockwise
        % which means they belong to the interior
        Iint = Iint + intOverPoly(f,polyCoords{k+4},Tol);
    end
end


end
