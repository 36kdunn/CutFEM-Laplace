function [extsln,intsln,x,y] = intExtFEM(N,G,flags)
%
% Main function for FEM solved separately on interior and exterior of
% interface

% Plot IC
% [uintX,uintY,uext,ugam,plot_ug,x,y] = packageOutput(N,u_prev);

gam1 = 1; % Dirichlet penalty enforcing parameter: u on boundary of Omega
gam2 = 1; % Velocity ghost penalty parameter
gam3 = 1; % Velocity 2nd order ghost penalty parameter
gam4 = 10; % Dirichlet penalty enforcing parameter: jump of u on Gamma

Tol = 1e-15;

% Define values from flags vector
intMeanZero  = flags(1);
extMeanZero  = flags(2);

% Store number of points given on Gamma
NGam = length(G); % NUMBER OF POINTS ON GAMMA PLUS 1 (FIRST COUNTED TWICE)

% Get connectivity table T and table of vertex coordinates V
[T,V] = connectivityTableQ2(N);

% Get list of curve segments and global indices of curve segments
[allCoords,segData] = findInterfaceIntersections(G,N);

% Find the elements that are intersected by the curve
T = colorSquares(T,N,allCoords,segData);

if N<32
    % warning('plotStuff is on!!!');          plotStuff(T,G,V,N)
end

% Get local indices of squares
[extLocalInds,intLocalInds,Next,Nint] = getLocalIndices(T,N);

% Get matrices to be added or used in Gaussian quadrature to avoid
% excessive/unnecessary revaluation of functions
mats = getMatsQ2(N);

% STRUCTURED WITH INTERIOR FIRST!!!
A = zeros(Nint + Next + 2);

% % % Initialize RHS vector
b = zeros(Nint + Next + 2,1);

% Define coordinates on [-1,1] for two-point Gaussian quadrature
gaussPts = [-1/sqrt(3), 1/sqrt(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Iterate over squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sqInd = 1:size(T,1)
    
    switch T(sqInd,10)
        
        case 1 % Exterior square
            % Get local indicies
            extInds =  Nint + extLocalInds(T(sqInd,1:9),1);
            
            x0 = V(T(sqInd,1),1);
            y0 = V(T(sqInd,1),2);
            
            % % % LHS contributions:
            % Contribution to exterior bulk from u - div(eps(u))
            A(extInds,extInds) = A(extInds,extInds) + mats.gradugradv(1:9,1:9);
            
            % % % RHS contributions:
            gaussX = 0.5*(gaussPts/N + 2*x0 + 1/N);
            gaussY = 0.5*(gaussPts/N + 2*y0 + 1/N);
            
            area = (1/N)^2;
            
            b(extInds,1) = b(extInds,1) + (area/4)*(fRHSext(gaussX(1,1),gaussY(1,1)).*mats.phiVals(:,1) +...
                fRHSext(gaussX(1,1),gaussY(1,2)).*mats.phiVals(:,2) + ...
                fRHSext(gaussX(1,2),gaussY(1,1)).*mats.phiVals(:,3) + ...
                fRHSext(gaussX(1,2),gaussY(1,2)).*mats.phiVals(:,4));
            
            IextOmegaDLHS = zeros(9);
            IextOmegaDRHS = zeros(9,1);
            
            if V(T(sqInd,4),2) == 0 % bottom edge of Omega
                
                n = [0,-1];
                
                % Contribution from Dirichlet enforcing penalty term
                inds = [1,4,7];
                IextOmegaDLHS(inds,inds) = ...
                    IextOmegaDLHS(inds,inds) + gam1*mats.uvOmega;
                
                % Contribution from unknown part of symmeterizing boundary term
                IextOmegaDLHS = IextOmegaDLHS - mats.graddotnvB - mats.graddotnvB';
                
                % Contribution from known part of Dirichlet penalty term
                IextOmegaDRHS = IextOmegaDRHS + gam1*gaussEdge(@(x,y) ...
                    fDomega(x,y).*basisFuncsVec(x,y,x0,y0,N),x0,x0+1/N,0,0)*N;
                
                % Contribution from known part of symmetrizing boundary term
                IextOmegaDRHS = IextOmegaDRHS - gaussEdge(@(x,y) ...
                    fgraddotn(x,0,x0,y0,N,n),x0,x0+1/N,0,0);
            end
            
            if V(T(sqInd,6),2) == 1 % top edge of Omega
                
                n = [0,1];
                
                % Contribution from Dirichlet enforcing penalty term
                inds = [3,6,9];
                IextOmegaDLHS(inds,inds) = ...
                    IextOmegaDLHS(inds,inds) + gam1*mats.uvOmega;
                
                % Contribution from unknown part of symmeterizing boundary term
                IextOmegaDLHS = IextOmegaDLHS - mats.graddotnvT - mats.graddotnvT';
                
                % Contribution from known part of Dirichlet penalty term
                IextOmegaDRHS = IextOmegaDRHS + gam1*gaussEdge(@(x,y) ...
                    fDomega(x,1).*basisFuncsVec(x,1,x0,y0,N),x0,x0+1/N,1,1)*N;
                
                % Contribution from known part of symmetrizing boundary term
                IextOmegaDRHS = IextOmegaDRHS - gaussEdge(@(x,y) ...
                    fgraddotn(x,1,x0,y0,N,n),x0,x0+1/N,1,1);
            end
            
            if V(T(sqInd,2),1) == 0 % left edge of Omega
                
                n = [-1,0];
                
                % Contribution from Dirichlet enforcing penalty term
                inds = [1,2,3];
                IextOmegaDLHS(inds,inds) = ...
                    IextOmegaDLHS(inds,inds) + gam1*mats.uvOmega;
                
                % Contribution from unknown part of symmeterizing boundary term
                IextOmegaDLHS = IextOmegaDLHS - mats.graddotnvL - mats.graddotnvL';
                
                % Contribution from known part of Dirichlet penalty term
                IextOmegaDRHS = IextOmegaDRHS + gam1*gaussEdge(@(x,y) ...
                    fDomega(0,y).*basisFuncsVec(0,y,x0,y0,N),0,0,y0,y0+1/N)*N;
                
                % Contribution from known part of symmetrizing boundary term
                IextOmegaDRHS = IextOmegaDRHS - gaussEdge(@(x,y) ...
                    fgraddotn(0,y,x0,y0,N,n),0,0,y0,y0+1/N);
            end
            
            if V(T(sqInd,8),1) == 1 % right edge of Omega
                
                n = [1,0];
                % Contribution from Dirichlet enforcing penalty term
                inds = [7,8,9];
                IextOmegaDLHS(inds,inds) = ...
                    IextOmegaDLHS(inds,inds) + gam1*mats.uvOmega;
                
                % Contribution from unknown part of symmeterizing boundary term
                IextOmegaDLHS = IextOmegaDLHS - mats.graddotnvR - mats.graddotnvR';
                
                % Contribution from known part of Dirichlet penalty term
                IextOmegaDRHS = IextOmegaDRHS + gam1*gaussEdge(@(x,y) ...
                    fDomega(1,y).*basisFuncsVec(1,y,x0,y0,N),1,1,y0,y0+1/N)*N;
                
                % Contribution from known part of symmetrizing boundary term
                IextOmegaDRHS = IextOmegaDRHS - gaussEdge(@(x,y) ...
                    fgraddotn(1,y,x0,y0,N,n),1,1,y0,y0+1/N);
                
            end
            
            A(extInds,extInds) = A(extInds,extInds) + IextOmegaDLHS;
            
            b(extInds,1) = b(extInds,1) + IextOmegaDRHS;
            
        case -1 % Interior square
            % Get local interior indices
            intInds = intLocalInds(T(sqInd,1:9),1);
            
            % % % LHS contributions:
            % Contribution to interior bulk from u - div(eps(u))
            A(intInds,intInds) = A(intInds,intInds) + mats.gradugradv;
            
            % % % RHS contributions:
            gaussX = 0.5*(gaussPts/N + 2*V(T(sqInd,1),1) + 1/N);
            gaussY = 0.5*(gaussPts/N + 2*V(T(sqInd,1),2) + 1/N);
            
            area = (1/N)^2;
            IRHSint = (area/4)*(fRHSint(gaussX(1,1),gaussY(1,1)).*mats.phiVals(:,1) +...
                fRHSint(gaussX(1,1),gaussY(1,2)).*mats.phiVals(:,2) + ...
                fRHSint(gaussX(1,2),gaussY(1,1)).*mats.phiVals(:,3) + ...
                fRHSint(gaussX(1,2),gaussY(1,2)).*mats.phiVals(:,4));
            
            b(intInds,1) = b(intInds,1) + IRHSint;
            
        case 0 % Interface square
            
            segs = cell(sum(segData(:,3) == sqInd),1);
            ind = 1;
            
            % Get local exterior indicies
            extInds = Nint + extLocalInds(T(sqInd,1:9),1);
            
            % Get local interior indices
            intInds = intLocalInds(T(sqInd,1:9),1);
            
            % Collect all curve segments in one cell array to pass to findPolygons
            for j = 1:length(segData)
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
            
            IextRHSbulk = zeros(9,1);
            IintRHSbulk = zeros(9,1);
            
            IextLHSbulk = zeros(9);
            IintLHSbulk = zeros(9);
            
            % Contribution to interface bulk from Laplacian
            for k = 1:4
                if ~isempty(polyCoords{k})
                    % First four cells in polyCoords are ordered clockwise
                    % which means they belong to the exterior
                    verts = zeros(size(polyCoords{k}));
                    [Nverts,~] = size(verts);
                    for l = 1:Nverts
                        verts(l,:) = polyCoords{k}(Nverts-l+1,:);
                    end
                    IextLHSbulk = IextLHSbulk + intOverPoly(@(x,y) LHSBulkExt(x,y,x0,y0,N),verts,Tol);
                    IextRHSbulk = IextRHSbulk + intOverPoly(@(x,y) RHSBulkExt(x,y,x0,y0,N),verts,Tol);
                end
                
                if ~isempty(polyCoords{k+4})
                    % Last four cells in polyCoords are ordered counterclockwise
                    % which means they belong to the interior
                    %                     polyCoords{5}=[x0, y0; x0+1/N, y0;x0+1/N,y0+1/N;x0, y0+1/N];
                    IintLHSbulk = IintLHSbulk + intOverPoly(@(x,y) LHSBulkInt(x,y,x0,y0,N),polyCoords{k+4},Tol);
                    IintRHSbulk = IintRHSbulk + intOverPoly(@(x,y) RHSBulkInt(x,y,x0,y0,N),polyCoords{k+4},Tol);
                end
            end
            
            % Add contributions to matrix & vector
            A(extInds,extInds) = A(extInds,extInds) + IextLHSbulk;
            
            b(extInds,1) = b(extInds,1) + IextRHSbulk;
            
            A(intInds,intInds) = A(intInds,intInds) + IintLHSbulk;
            
            b(intInds,1) = b(intInds,1) + IintRHSbulk;
            
    end % end switch
    
end % end loop over squares

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Iterate over segments of Gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if norm(G(1,:)-allCoords(1,:))==0
    Gind = 1;
else
    Gind = NGam-1;
end

for i = 1:length(segData)
    
    globalInds = T(segData(i,3),1:4);
    
    intInds = intLocalInds(T(segData(i,3),1:9),1);
    extInds = Nint + extLocalInds(T(segData(i,3),1:9),1);
    
    % Define coordinates of bottom-left corner of current square
    x0 = V(globalInds(1,1),1);
    y0 = V(globalInds(1,1),2);
    
    for j = segData(i,1):(segData(i,1)+segData(i,2)-2)
        
        % Compute the UNIT outward normal (outward to the exterior of Gamma)
        n = -[(allCoords(j+1,2)-allCoords(j,2)), -(allCoords(j+1,1)-allCoords(j,1))];
        n = n/norm(n);
        
        % % % Symmetrized jump terms on LHS
        % % Strain rate tensor term
        % LHS block contribution to jump condition; ex: +/-(eps(u),v)/2 +/- (eps(v),u)/2
        ILHSepsjump = zeros(9);
        ILHSepsjump = ILHSepsjump + gaussEdge(@(x,y) graddotnv(x,y,x0,y0,N,n),...
            allCoords(j,1),allCoords(j+1,1),allCoords(j,2),allCoords(j+1,2));
        
        A(intInds,intInds) = A(intInds,intInds) + 0.5*ILHSepsjump + 0.5*ILHSepsjump';
        A(intInds,extInds) = A(intInds,extInds) + 0.5*ILHSepsjump - 0.5*ILHSepsjump';
        A(extInds,intInds) = A(extInds,intInds) - 0.5*ILHSepsjump + 0.5*ILHSepsjump';
        A(extInds,extInds) = A(extInds,extInds) - 0.5*ILHSepsjump - 0.5*ILHSepsjump';
        
        % % % Penalty term enforcing jump of u
        %  LHS term: penalty term on jump of u (gam6/h)*([u],[v])_Gamma
        ILHSujumppenalty = zeros(9);
        ILHSujumppenalty = ILHSujumppenalty + gam4*N*gaussEdge(@(x,y) ...
            basisFuncsVec(x,y,x0,y0,N)*basisFuncsVec(x,y,x0,y0,N)',...
            allCoords(j,1),allCoords(j+1,1),allCoords(j,2),allCoords(j+1,2));
        
        A(extInds,extInds) = A(extInds,extInds) + ILHSujumppenalty;
        A(intInds,extInds) = A(intInds,extInds) - ILHSujumppenalty;
        A(extInds,intInds) = A(extInds,intInds) - ILHSujumppenalty;
        A(intInds,intInds) = A(intInds,intInds) + ILHSujumppenalty;
        
        % % % RHS term: penalty term on jump of u
        % (gam6/h)*(fjumpuX/fjumpuY,[v])_Gamma
        IRHSujumppenalty = zeros(9,1);
        IRHSujumppenalty = IRHSujumppenalty + gam4*N*gaussEdge(@(x,y) ...
            fjumpu(x,y)*basisFuncsVec(x,y,x0,y0,N),...
            allCoords(j,1),allCoords(j+1,1),allCoords(j,2),allCoords(j+1,2));
        
        % % % RHS term: jump of u (result of symmetrization of LHS)
        % (fjumpuX/fjumpuY,{eps(v)n})
        % Subtract 0.5*IRHSujump from both interior and exterior RHS
        IRHSujump = zeros(9,1);
        IRHSujump = IRHSujump + gaussEdge(@(x,y) RHSGamJumpU(x,y,x0,y0,N,n),...
            allCoords(j,1),allCoords(j+1,1),allCoords(j,2),allCoords(j+1,2));
        
        
        % % % RHS term: jump of (grad(u)-p)n (result of symmetrization of LHS)
        % (fjumpgraduX/fjumpgraduY,{v})
        % Add 0.5*IRHSjumpgradu to both interior and exterior
        IRHSjumpgradu = zeros(9,1);
        IRHSjumpgradu = IRHSjumpgradu + gaussEdge(@(x,y) RHSGamJumpGradU(x,y,x0,y0,N),...
            allCoords(j,1),allCoords(j+1,1),allCoords(j,2),allCoords(j+1,2));
        
        % Add contributions of Dirichlet and Neumann data on Gamma to RHS vectors
        b(extInds,1) = b(extInds,1) - 0.5*IRHSujump + 0.5*IRHSjumpgradu + IRHSujumppenalty;
        b(intInds,1) = b(intInds,1) - 0.5*IRHSujump + 0.5*IRHSjumpgradu - IRHSujumppenalty;
        
        % Update Gind if necessary
        if norm(allCoords(j+1,:) - G(Gind+1,:))==0
            if Gind < NGam-1
                Gind = Gind+1;
            else
                Gind = 1;
            end
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Ghost penalty term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add contributions for ghost penalty term
A = ghostPenalty(T,N,A,0,intLocalInds,gam2,gam3,mats,-1); % ghost penalty for interior
A = ghostPenalty(T,N,A,Nint,extLocalInds,gam2,gam3,mats,1); % ghost penalty for exterior

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Add Lagrange multiplier (if necessary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If necessary, force impose the mean-zero condition on the interior/exterior solution
% if intMeanZero+extMeanZero>0
[Wext,Wint] = getWeights(T,V,N,allCoords,segData,extLocalInds,intLocalInds,Tol);

% Compute row/column index of first weight vector
Windex1 = Nint+Next+1;

if intMeanZero == 1
    
    % Weights for interior
    A(1:Nint,Windex1) = Wint;
    A(Windex1,1:Nint) = Wint';
    
else
    A(Windex1,Windex1) = 1;
end

if extMeanZero == 1
    % Weights for exterior
    A(Nint + (1:Next),Windex1+1) = Wext;
    A(Windex1+1,Nint + (1:Next)) = Wext';
else
    A(Windex1+1,Windex1+1) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Solve the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b(1:2*Nint+Next+NGam-1,1) = b(1:2*Nint+Next+NGam-1,1) + Ivuprev;
u = A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Error calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uafunc.ext = @(x,y) x.*exp(-x.*y);
uafunc.int = @(x,y) x.*exp(-x.*y);

uafunc.extXx = @(x,y) exp(-x.*y) - x.*y.*exp(-x.*y);
uafunc.extXy = @(x,y) -(x.^2).*exp(-x.*y);
uafunc.intXx = @(x,y) exp(-x.*y) - x.*y.*exp(-x.*y);
uafunc.intXy = @(x,y) -(x.^2).*exp(-x.*y);
%{
% Compute L2 error of solutions on respective domains
[EintuXL2,EintuYL2,EextuXL2,EextuYL2,EintpL2,EextpL2] = compL2Error(u,uafunc,T,Tp,V,N,Nint,Nintp,Next,...
    allCoords,segData,extLocalInds,intLocalInds,extLocalIndsP,intLocalIndsP,Tol);

% Compute H1 error of solutions on respective domains
[EintuXH1,EintuYH1,EextuXH1,EextuYH1,EintpH1,EextpH1] = compH1Error(u,uafunc,T,Tp,V,N,Nint,Nintp,Next,...
    allCoords,segData,extLocalInds,intLocalInds,extLocalIndsP,intLocalIndsP,Tol);

% Compute Linfty error of solutions on respective domains
[EintuX,EintuY,EextuX,EextuY,Eintp,Eextp,xyExt,xyInt] = compLinftyError(u,uafunc,T,Tp,V,Vp,N,Nint,Nintp,Next,...
    allCoords,segData,extLocalInds,intLocalInds,extLocalIndsP,intLocalIndsP,Tol);

% Compute W1,infty error of solutions on respective domains
[EintuXW1,EintuYW1,EextuXW1,EextuYW1,EintpW1,EextpW1,xyExtW1,xyIntW1] = compW1inftyError(u,uafunc,T,Tp,V,Vp,N,Nint,Nintp,Next,...
    allCoords,segData,extLocalInds,intLocalInds,extLocalIndsP,intLocalIndsP,Tol);

format shorte


disp('                 L2 Error')
disp('     u_x     |    u_y    |    p    ')
disp('------------------------------------')
% disp([sqrt(EintuXL2^2+EextuXL2^2),sqrt(EintuYL2^2+EextuYL2^2),sqrt(EintpL2^2 + EextpL2^2)])
disp([sqrt(EintuXL2^2 + EextuXL2^2 + EintuYL2^2 + EextuYL2^2),sqrt(EintpL2^2 + EextpL2^2)])

disp('                 H1 Error')
disp('     u_x     |    u_y    |    p    ')
disp('------------------------------------')
% disp([sqrt(EintuXH1^2+EextuXH1^2),sqrt(EintuYH1^2+EextuYH1^2),sqrt(EintpH1^2 + EextpH1^2)])
disp([sqrt(EintuXH1^2 + EintuYH1^2 ) + sqrt(EextuXH1^2 + EextuYH1^2),sqrt(EintpH1^2 + EextpH1^2)])

disp('              Linfty Error')
disp('     u_x     |    u_y    |    p    ')
disp('------------------------------------')
% disp([max(EintuX,EextuX),max(EintuY,EextuY),max(Eintp,Eextp)])
disp([max(EintuX+EintuY,EextuX+EextuY),max(Eintp,Eextp)])


disp('            W1infty Error')
disp('     u_x     |    u_y    |    p    ')
disp('------------------------------------')
% disp([EintuXW1,EextuXW1,EintuYW1,EextuYW1,EintpW1,EextpW1])
% disp([max(EintuXW1,EextuXW1),max(EintuYW1,EextuYW1),max(EintpW1,EextpW1)])
disp([max(EintuXW1+EintuYW1,EextuXW1+EextuYW1),max(EintpW1,EextpW1)])




format short
%}

ua = zeros(size(b));
uatempext = zeros(Next,1);
uatempint = zeros(Nint,1);

for i = 1:(2*N+1)^2
    x = V(i,1);
    y = V(i,2);
    if extLocalInds(i,1)>0
        % Linear Example
        %         uatempxext(extLocalInds(i,1),1) = 3+x+y-4;
        %         uatempyext(extLocalInds(i,1),1) = x+2*y-1.5;
        
        % Exponential/Trig Example
        uatempext(extLocalInds(i,1),1) = uafunc.ext(x,y);
        
        % RBM Test Example
        %         uatempxext(extLocalInds(i,1),1) = 2-2*y;
        %         uatempyext(extLocalInds(i,1),1) = 1+2*x;
        
    end
    
    if intLocalInds(i,1)>0
        % Linear Example
        %         uatempxint(intLocalInds(i,1),1) = 3+x+y-4;
        %         uatempyint(intLocalInds(i,1),1) = x+2*y-1.5;
        
        % Exponential/Trig Example
        uatempint(intLocalInds(i,1),1) = uafunc.int(x,y);
        
        % RBM Test Example
        %         uatempxint(intLocalInds(i,1),1) = 2-2*y;
        %         uatempyint(intLocalInds(i,1),1) = 1+2*x;
        
    end
    
end

%Interior actual solution
ua(1:Nint,1) = uatempint;

% Exterior actual solution
ua(Nint + (1:Next),1) = uatempext;

u=abs(u-ua);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Package output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = kron(ones(2*N+1,1),0:1/(2*N):1);
x = y';

extsln.u = zeros(2*N+1);
intsln.u = zeros(2*N+1);

for i = 1:(2*N+1)^2
    xloc = mod(i,2*N+1);
    if xloc == 0
        xloc = 2*N+1;
    end
    yloc = (i-xloc)/(2*N+1)+1;
    
    if intLocalInds(i,1) > 0
        intsln.u(xloc,yloc) = u(intLocalInds(i,1),1);
    else
        intsln.u(xloc,yloc) = NaN;
    end
    
    if extLocalInds(i,1) > 0
        extsln.u(xloc,yloc) = u(Nint + extLocalInds(i,1),1);
    else
        extsln.u(xloc,yloc) = NaN;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Display plots (if uncommented)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% {
figure();
subplot(1,2,1)
surf(x,y,extsln.u)
subplot(1,2,2)
surf(x,y,intsln.u)
%}
end

