function [x,w,P] = gll_nodes(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% gll_nodes.m
%
% Computes the Gauss-Lobatto-Legendre (GLL) nodes, weights and the Vandermonde 
% matrix. The GLL nodes are the zeros of (1-x^2)*P'_N(x), where P'_N(x) is the
% first derivative of the Legendre polynomial of order N. The derived set of
% points and weights is useful for numerical integration and spectral methods. 
%
% Reference on GLL nodes and weights: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of points
N1 = N + 1;

% Use the Gauss-Lobatto-Chebyshev nodes as an initial guess
x = cos(pi*(0:N)/N)';

% Legendre Vandermonde Matrix
P = zeros(N1,N1);

% Compute P_(N) using the recursion relation.
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold = 2*ones(size(x));

while max(abs(x-xold)) > eps

    xold = x;
        
    P(:,1) = 1;
	P(:,2) = x;
    
    for k = 2:N
        P(:,k+1) = ((2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1))/k;
    end
     
    x = xold - (x.*P(:,N1)-P(:,N) )./( N1*P(:,N1));
             
end

w = 2./(N*N1*P(:,N1).^2);
x = -x;

end