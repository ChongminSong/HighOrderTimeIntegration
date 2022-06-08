function [p,q,r,C] = TimeIntgCoeff(M,w)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Chongmin Song & Sascha Eisentraeger (13/10/2021)
%
% SPDX-License-Identifier: MIT
%
%-------------------------------------------------------------------------%
% Input:
% order:  Order of the Pade expansion (11,22,33,44,55)
% w:      Weight factor for the mixed Pade expansion (spectral radius)
%
% Output:
% q:      Coefficients of the numerator polynomial of the Pade expansion
% p:      Coefficients of the denominator polynomial of the Pade expansion
% r:      Roots of the denominator polynomial of the Pade expansion
% C:      Multipliers for the force vector (C_i-matrices)
%-------------------------------------------------------------------------%

M = M/11;

%% Diagonal Pade expansion
L = M;
% Compute the coefficients of the numerator polynomial
ii = 0:L; 
p1 = (factorial(M+L-ii))./(factorial(ii).*factorial(L-ii));
% Compute the coefficients of the denominator polynomial
ii = 0:M; 
q1 = (factorial(M+L-ii)).*((-1).^ii)./(factorial(ii).*factorial(M-ii))*...
      factorial(M)/factorial(L);

%% Mixed order Pade expansion
L = M-1;
% Compute the coefficients of the numerator polynomial
ii = 0:L; 
p2 = (factorial(M+L-ii))./(factorial(ii).*factorial(L-ii));
% Compute the coefficients of the denominator polynomial
ii = 0:M; 
q2 = (factorial(M+L-ii)).*((-1).^ii)./(factorial(ii).*factorial(M-ii))*...
      factorial(M)/factorial(L);
  
%% Weighted Pade expansion
p = w*p1 + (1-w)*[p2 zeros(M-L)];
q = w*q1 + (1-w)*q2;
% Compute the roots of the denominator polynomial
r = roots(q(end:-1:1));
r(imag(r)<-1.d-6) = [];
[~, idx] = sort(imag(r));
r = r(idx);

% Compute coefficients for the C_i-matrices
tmp = p - q;
C = zeros(M+1,M);
C(1,:) = tmp(2:end);
f1 = 1;
f2 = 1;
for ii = 1:M
    f1 = -f1/2;
    f2 = -f2;
    tmp = ii*[C(ii,:) zeros(M-L)] + f1*(p-f2*q);
    C(ii+1,:) = tmp(2:end);
end

end