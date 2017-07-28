function [extLocalInds, intLocalInds,Next,Nint] = getLocalIndices(T,N)
%
% Obtain local index numbering for interior and exterior nodes. These
% tables are necessary for the indexing of the matrix and recovering the
% solution after matrix solve. The tables are vectors of length (2N+1)^2
% which have a nonzero entry containing the local index if the node at the 
% global index is not contained in the corresponding mesh and contains a 
% zero index otherwise.
% Inputs:
%   T = Connectivity table
%   N = Number of squares along one edge of Omega.
% Outputs:
%   extLocalInds = vector/table of local indexing for exterior nodes
%   intLocalInds = vector/table of local indexing for interior nodes
%   Next = number of exterior DOF
%   Nint = number of interior DOF

% Initialize output vectors and indices
intLocalInds = zeros((2*N+1)^2,1);
extLocalInds = intLocalInds;
intInd = 1;
extInd = 1;

for i = 1:N^2 % length(T)
    
    if T(i,10) >= 0
        % If node is in exterior mesh
        for j = 1:9
            % Add all indices to local index table
            if extLocalInds(T(i,j),1) == 0
                extLocalInds(T(i,j),1) = extInd;
                extInd = extInd + 1;
            end
        end
    end
    
    if T(i,10) <= 0
        % If node is in interior mesh
        for j = 1:9
            % Add all indices to local index table
            if intLocalInds(T(i,j),1) == 0
                intLocalInds(T(i,j),1) = intInd;
                intInd = intInd + 1;
            end
        end
    end
end
% Subtract 1 to obtain the total number of degrees of freedom in each mesh
Next = extInd - 1;
Nint = intInd - 1;
end