function [dsp,vel,acc] = PadeSolver11(order,rho_inf,ft,K,M,C,F,U0,V0,...
                         saveDOF,nnzIC,flag_wb)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Chongmin Song & Sascha Eisentraeger (13/10/2021)
%
% SPDX-License-Identifier: MIT
%
%-------------------------------------------------------------------------%
% Input
% order:       	11
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

[pcoe,qcoe,rs,cf] = TimeIntgCoeff(order,rho_inf);
r = rs(1);
Kd = sparse((r*r)*M + r*C + K);
dKd = decomposition(Kd);

nDOF = size(K,1);

f  = [Fb; zeros(nDOF,1)];
Pf = f*cf';

npt = size(ft,2);
Pf = Pf(:,1:npt);
pt = 0.5.^(0:npt-1)';

if nnzIC ~= 0
	x = [V0; U0];
    Ax1 = SolverPadeAx(dM,K,C,x);
	
    % Initial conditions
    dsp(:,1) = U0(saveDOF);
    vel(:,1) = V0(saveDOF);
	
	% Initial acceleration --> normalized with dt^2
	fn = ft(1,:)*((-0.5).^(0:npt-1))'; % force at t=0
    acc(:,1) = Ax1(saveDOF) + fn*Fb(saveDOF);
else
	x = zeros(2*nDOF,1);
    Ax1 = x;
end

if rho_inf == 1
    strg = '2nd-order Pade-based time integrator (implicit): ';
else
    strg = '1st-order Pade-based time integrator (implicit): ';
end

if flag_wb
    % Waitbar
    h = waitbar(0,strg);

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ndt/1000);
end

disp('Solve the equations of motion using the Pade-based time integrator');
disp('');
disp('---------------------------------------------------------------------');
disp('');

for it = 2:ndt
	is = it - 1;
	
    b = [x Ax1]*pcoe' + Pf*ft(is,:)';
    ctmp = r*M*b(1:nDOF) - K*b(nDOF+1:2*nDOF);
    tmp = dKd\ctmp;
    x = [tmp; (tmp+b(nDOF+1:2*nDOF))/r];

    Ax1 = (b - qcoe(1)*x)/qcoe(2);
    
    % store responses for output
    vel(:,it) = x(saveDOF);
    dsp(:,it) = x(nDOF+saveDOF);
	fn = ft(is,:)*pt;
    acc(:,it) = Ax1(saveDOF) + fn*Fb(saveDOF);
    
    if flag_wb
        % Waitbar
        if mod(it,nUpdate) == 0
            waitbar(it/ndt,h,[strg num2str(it/ndt*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end