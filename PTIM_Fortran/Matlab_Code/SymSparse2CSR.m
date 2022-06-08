function [ia, ja, aK, aM, aC] = SymSparse2CSR(K, M, C)
% Convert MATLAB storage of a symmetric sparse matrix to 
% CSR (Compressed Sparse Row) format

% Number of degrees of freedom
n = size(K,1);

% Different output depending on the number of input arguments
switch nargin
    case (1)
        [ja, ia, aK] = find(tril(K));
        aM = []; aC = [];
    case (2)
        [ja, ia] = find(abs(tril(K))+abs(tril(M)));
        id = (ia-1)*n+ja;
        aK = full(K(id)); aM = full(M(id)); aC = [];
    case (3)
        [ja, ia] = find(abs(tril(K))+abs(tril(M))+abs(tril(C)));
        id = (ia-1)*n+ja;
        aK = full(K(id)); aM = full(M(id)); aC = full(C(id));
end
ia = cumsum([ 1; accumarray(ia,1,[n 1])]);

end