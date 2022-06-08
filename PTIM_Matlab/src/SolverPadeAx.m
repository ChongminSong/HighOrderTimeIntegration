function y = SolverPadeAx(dM,K,C,x)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Chongmin Song & Sascha Eisentraeger (13/10/2021)
%
% SPDX-License-Identifier: MIT
%
%-------------------------------------------------------------------------%
% Input:
% dM:	factorized mass matrix, decomposition object
% K:  	stiffness matrix
% C:  	damping matrix
% x:   	state-space vector
%
% Output:
% y:    Result vector
%-------------------------------------------------------------------------%%

% Number of degrees of freedom (x is a vector in state-space)
n = length(x)/2;

% Auxiliary variable
x1 = x(1:n);
x2 = x(n+1:2*n);
tmp = C*x1 + K*x2;

% Result
y = [-dM\tmp ; x(1:n)];

end