function [dsp,vel,acc] = PadeSolver33(order,rho_inf,ft,K,M,C,F,U0,V0,...
                         saveDOF,nnzIC,flag_wb)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Chongmin Song & Sascha Eisentraeger (13/10/2021)
%
% SPDX-License-Identifier: MIT
%
%-------------------------------------------------------------------------%
% Input
% order:        33
% rho_inf:      spectral radius of the amplification matrix
% ft:           external loading; force matrix for the Pade time
%               integration scheme
% K:            constrained stiffness matrix
% M:            constrained mass matrix
% C:            constrained damping matrix
% F:            constrained force vector
% U0:           initial displacements
% V0:           initial velocities
% saveDOF:    	the DOF for output response
% nnzIC:        number of DOFs that have non-zero initial conditions
% flag_wb:      (de-)active the waitbar
%
% Output
% dsp:          displacement vector
% vel:          velocity vector
% acc:          acceleration vector
%-------------------------------------------------------------------------%
%%

% number of time steps
ndt = size(ft,1) + 1;

% initializing variables storing response history for output
dsp = zeros(length(saveDOF),ndt); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

dM = decomposition(M);
Fb = (dM\F);

% integration coefficients
[pcoe,qcoe,rs,cf] = TimeIntgCoeff(order,rho_inf);

% matrix decomposition 1
r1 = real(rs(1));
Kd1 = sparse((r1*r1)*M + r1*C + K);
dKd1 = decomposition(Kd1);

% matrix decomposition 2
r2 = rs(2);
Kd2 = sparse((r2*r2)*M + r2*C + K);
[L,U,p,q] = lu(Kd2,'vector');

nDOF = size(K,1);

f  = [Fb; zeros(nDOF,1)];
Af1 = SolverPadeAx(dM,K,C,f);
Af2 = SolverPadeAx(dM,K,C,Af1);
Pf = [f Af1 Af2]*cf';

npt = size(ft,2);
Pf = Pf(:,1:npt);
pt = 0.5.^(0:npt-1)';

if nnzIC ~= 0
	x = [V0; U0];
    Ax1 = SolverPadeAx(dM,K,C,x);
    Ax2 = SolverPadeAx(dM,K,C,Ax1);
    Ax3 = SolverPadeAx(dM,K,C,Ax2);
	
    % Initial conditions
    dsp(:,1) = U0(saveDOF);
    vel(:,1) = V0(saveDOF);
	
	% Initial acceleration --> normalized with dt^2
	fn = ft(1,:)*((-0.5).^(0:npt-1))'; % force at t=0
    acc(:,1) = Ax1(saveDOF) + fn*Fb(saveDOF);
else
	x = zeros(2*nDOF,1);
    Ax1 = x;
    Ax2 = x;
    Ax3 = x;
end

if rho_inf == 1
    strg = '6th-order Pade-based time integrator (implicit): ';
else
    strg = '5th-order Pade-based time integrator (implicit): ';
end

if flag_wb
    % Waitbar
    h = waitbar(0,strg);

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ndt/1000);
end

for it = 2:ndt
	is = it - 1;
    
	b = [x Ax1 Ax2 Ax3]*pcoe' + Pf*ft(is,:)';
    tmp = dKd1\(M*(r1*b(1:nDOF)) - K*b(nDOF+1:2*nDOF));
    x1 = [tmp; (tmp+b(nDOF+1:2*nDOF))/r1];
    
    ctmp = r2*M*x1(1:nDOF) - K*x1(nDOF+1:2*nDOF);
    tmp = U\(L\ctmp(q));
    tmp(p) = tmp;
    y = [tmp; (tmp+x1(nDOF+1:2*nDOF))/r2];
    x = -imag(y)/imag(r2);
    
    Ax1 = -imag(r2*y)/imag(r2);
    Ax2 = x1 - real(r2*conj(r2))*x + 2*real(r2)*Ax1;
    Ax3 = (b - [x Ax1 Ax2]*qcoe(1:3)')/qcoe(4);
    
    % store responses for output
    vel(:,it) = x(saveDOF);
    dsp(:,it) = x(nDOF+saveDOF);
	fn = ft(is,:)*pt;
    acc(:,it) = Ax1(saveDOF) + fn*Fb(saveDOF);
    
    if flag_wb
        % Waitbar
        if mod(it,nUpdate) == 0
            waitbar((it-1)/(ndt-1),h,[strg num2str((it-1)/(ndt-1)*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end