function [t,dsp,vel,acc] = PTIM(order,rho_inf,ndt,dt,ft,K,M,C,F,U0,V0,...
                           saveDOF,flag_wb)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Chongmin Song & Sascha Eisentraeger (13/10/2021)
%
% SPDX-License-Identifier: MIT
%
%-------------------------------------------------------------------------%
% Input
% order:        order of Pade expansion: = 11, 22, 33, 44, 55
% rho_inf:      Spectral radius of the amplification matrix
% ndt:       	number of time steps
% dt:           time steps size
% ft:           external loading; force matrix for the Pade time
%               integration scheme
% K:            constrained stiffness matrix
% M:            constrained mass matrix
% C:            constrained damping matrix
% F:            constrained force vector
% U0:           initial displacements
% V0:           initial velocities
% saveDOF:    	the DOF for output response
% flag_wb:      (de-)active the waitbar
%
% Output:
% t:            time vector
% dsp:          displacement vector
% vel:          velocity vector
% acc:          acceleration vector
%-------------------------------------------------------------------------%
%% Initial conditions
nDOF = size(K,1);  % number of DOFs
if isempty(U0)
    U0 = zeros(nDOF,1);
    nnzU0 = 0;
else
    % Number non-zero terms
    nnzU0 = nnz(U0);
end
if isempty(V0)
    V0 = zeros(nDOF,1);
    nnzV0 = 0;
else
    % Number non-zero terms
    nnzV0 = nnz(V0);
end
% Number of non-zero terms in the initial conditions
nnzIC = nnzU0 + nnzV0;

%% exclude fixed DOFs in the analysis
dt2 = dt*dt;
% Stiffness, damping, and mass matrices
K = dt2*K;
C = dt*C;
% Force vector
F = dt2*F;
% Initial velocity
V0 = dt*V0;

% Command window output: Time step size
fprintf('Time step size used: %6.3e\n',dt);

switch order
    case 11
        [dsp,vel,acc] = PadeSolver11(order,rho_inf,ft,K,M,C,F,U0,V0,...
                        saveDOF,nnzIC,flag_wb);
	case 22
        [dsp,vel,acc] = PadeSolver22(order,rho_inf,ft,K,M,C,F,U0,V0,...
                        saveDOF,nnzIC,flag_wb);
    case 33
        [dsp,vel,acc] = PadeSolver33(order,rho_inf,ft,K,M,C,F,U0,V0,...
                        saveDOF,nnzIC,flag_wb);
    case 44
        [dsp,vel,acc] = PadeSolver44(order,rho_inf,ft,K,M,C,F,U0,V0,...
                        saveDOF,nnzIC,flag_wb);
    case 55
        [dsp,vel,acc] = PadeSolver55(order,rho_inf,ft,K,M,C,F,U0,V0,...
                        saveDOF,nnzIC,flag_wb);
    otherwise
        % The variable order is incorrectly defined
        error("The value for order is not supported!");
end

dsp = dsp';
vel = vel'/dt;
acc = acc'/(dt*dt);
t = (0:ndt)'*dt;

end