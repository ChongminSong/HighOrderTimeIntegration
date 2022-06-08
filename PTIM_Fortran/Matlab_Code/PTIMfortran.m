function [t,dsp,vel,acc] = PTIMfortran(pathExe,order,rho_inf,dt,ft,K,M,...
                                       C,F,U0,V0,saveDOF)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Chongmin Song & Sascha Eisentraeger (01/06/2022)
%
% SPDX-License-Identifier: MIT
%
%-------------------------------------------------------------------------%
% Input
% pathExe:      Complete path to the executable Fortran code
% order:        order of Pade expansion: = 11, 22, 33, 44, 55
% rho_inf:      Spectral radius of the amplification matrix
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

% Create input for the FORTRAN solver
inputPath = 'c:\tmp\';
% Create temporary folder if it does not exist
if ~exist(inputPath, 'dir')
	mkdir(inputPath)
end

% Control file
fname = [inputPath 'cntrl.txt'];
fileID = fopen(fname,'w');
fprintf(fileID,'%s\n', "porder rho_infty");
fprintf(fileID,'%d  %g\n',order,rho_inf);
fnameKM = 'KM.txt';
fprintf(fileID,'KM %s\n',fnameKM);
fnameFR = 'Force.txt';
fprintf(fileID,'FR %s\n',fnameFR);
if nnzIC > 0
	fnameIC = 'InitCond.txt';
	fprintf(fileID,'IC %s\n',fnameIC);
end
fnameTH = 'TimeHistrory.txt';
fprintf(fileID,'TH %s\n',fnameTH);
fnameOC = 'OutputControl.txt';
fprintf(fileID,'OC %s\n',fnameOC);
fnameRS = 'response01.txt';
fprintf(fileID,'RS %s\n',fnameRS);
fprintf(fileID,'END\n');
fclose(fileID);

% Convert matrices to CSR format
[ia,ja,Kc,Mc,Cc] = SymSparse2CSR(K,M,C);
% Save matrices in an ASCII-file
fileID = fopen([inputPath fnameKM],'w');
OutputMatrixCSR(fileID,ia,ja,Kc,Mc,Cc);
fclose(fileID);

% Save force vector in an ASCII-file
fileID = fopen([inputPath fnameFR],'w');
fprintf(fileID,'%s\n', "F");
fprintf(fileID,'%-23.16g\n',full(F));
fclose(fileID);

% Save initial conditions in an ASCII-file
if nnzIC > 0
	fnameIC = 'InitCond.txt';
	fileID = fopen([inputPath fnameIC],'w');
	fprintf(fileID,'%s\n', "u0");
	fprintf(fileID,'%-23.16g\n',U0);
	fprintf(fileID,'%s\n', "v0");
	fprintf(fileID,'%-23.16g\n',V0);
end

% List of degrees of freedom that are saved
fileID = fopen([inputPath fnameOC],'w');
fprintf(fileID,'%s\n', "pDOF");
fprintf(fileID,'%d\n', length(saveDOF));
fprintf(fileID,'%-23.16g\n',saveDOF);
fclose(fileID);

% Save force vector in an ASCII-file
fileID = fopen([inputPath fnameTH],'w');
fprintf(fileID,'%s\n', "dt");
fprintf(fileID,'%-23.16g\n',dt);
fprintf(fileID,'%s\n', "dim of ft");
fprintf(fileID,'%d %d\n',size(ft));
fprintf(fileID,'%s\n', "ft");
fprintf(fileID,'%-23.16g\n',full(ft));
fclose(fileID);

% Call FORTRAN solver
command = ['"' pathExe '" c:\tmp\'];
t0 = tic;
system(command,'-echo');
toc(t0);

% Save results in an ASCII-file
outfname = [inputPath 'response01.txt'];
outfileID = fopen(outfname,'r');
fgets(outfileID);
[ns_f] = fscanf(outfileID, "%d\n",1);
fgets(outfileID);
[npDOF_f] = fscanf(outfileID, "%d\n",1);
fgets(outfileID);
t = fscanf(outfileID, "%g\n",ns_f);
fgets(outfileID);
dsp = fscanf(outfileID, "%g\n",ns_f*npDOF_f);
dsp = reshape(dsp,ns_f,[]);
fgets(outfileID);
vel = fscanf(outfileID, "%g\n",ns_f*npDOF_f);
vel = reshape(vel,ns_f,[]);
fgets(outfileID);
acc = fscanf(outfileID, "%g\n",ns_f*npDOF_f);
acc = reshape(acc,ns_f,[]);
fclose(outfileID);

% Save computational time in an ASCII-file
fname = [inputPath 'timing.txt'];
fileID = fopen(fname,'r');
[FortranTime] = fscanf(fileID, "%g",1);
disp(['Time (Fortran) for time integration = ', num2str(FortranTime), 's']);
fclose(fileID);

end