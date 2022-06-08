%% Example call for the high-order Pade-base time integration method
% Reference: Song, C. and Eisentraeger, S.: High-order implicit time 
% integration scheme based on Pade expansions, arXiv:2103.12282

%% Minimal Example: Three degree-of-freedom spring system
% Three masses (m1 to m3) are connected by two springs (k1 and k2)
% System: m1---k1---m2---k2---m3

%% 
clearvars; clc; close all

% Load all folders into the Matlab search path
path1 = genpath(['.' filesep 'Fortran_Code']);
path2 = genpath(['.' filesep 'Matlab_Code']);
addpath(path1,path2);

%% Path to executable Fortran code
pathExe = '.\Fortran_Code\x64\Release\HighOrderTimeIntg.exe';

%% Parameters
% Spring constant
k1 = 5;
k2 = 1;
% Mass
m1 = 1;
m2 = 2;
m3 = 1;
% Number of degrees of freedom
nDOF = 3;
% Initial conditions (displacement and velocity at t=0)
U0 = zeros(nDOF,1);
V0 = zeros(nDOF,1);
% Excitation
amp = 1;                            % amplitude
omega_ex = 1.2;                     % excitation frequency
signal = @(t) cos(omega_ex*t);      % excitation signal

%% Settings for the time integrator
tsim = 15;      % Simulation time
dt = 1e-2;      % Time step size/increment
rho_inf = 1.0;	% Spectral radius of the time integration scheme
order = 33;     % Order of the Pade expansion (11,22,33,44,55)

%% Analysis
% Mass matrix
M = diag([m1,m2,m3]);
% Damping matrix
C = zeros(nDOF);
% Stiffness matrix
K = [k1,-k1,0; -k1,k1+k2,-k2; 0,-k2,k2];
% Force vector
F = [0;1;0];    % external force acts on mass m3
% Fixed degrees of freedom
fixedDOF = 1;
% Active degrees of freedom
activeDOF = setdiff((1:nDOF)',fixedDOF);
% Save results for particular degrees of freedom
saveDOF = [2;3];

% Constrained system matrices and force vector
Mc = M(activeDOF,activeDOF);
Cc = C(activeDOF,activeDOF);
Kc = K(activeDOF,activeDOF);
Fc = F(activeDOF);
U0 = U0(activeDOF);
V0 = V0(activeDOF);

% Identify DOF for output in the reduced system
vecID = zeros(nDOF,1);
vecID(activeDOF) = 1:length(activeDOF);
saveDOF = vecID(saveDOF);
saveDOF(saveDOF == 0) = [];

% Number of time steps
ndt = ceil(round(tsim/dt,3,'decimal'));
% Adjust the time step size if necessary
dt = tsim/ndt;
% Check if an admissible value for the order has been defined
if mod(order,11) ~= 0
    order = floor(order/11)*11;
end
% Check if rho_inf is in the range [0,1]
if rho_inf > 1
    rho_inf = 1.0;
end
if rho_inf < 0
    rho_inf = 0.0;
end

%% Eigenvalue analysis
[MS,EV] = eigs(Kc,Mc);
f_max = max(sqrt(diag(EV)))/(2*pi);     % Largest natural frequency
T_min = 1/f_max;                        % Smallest natural period
nP = T_min/dt;                          % Time steps per smallest period

%% Time integration
% Pre-compute the amplitude of the force vector for a polynomial variation
% of the force within a time step
ft = ForcePTIM(order,signal,ndt,dt);

% Call the time integrator
[T,U,V,A] = PTIMfortran(pathExe,order,rho_inf,dt,ft,Kc,Mc,Cc,Fc,U0,V0,...
                        saveDOF);

%% Visualization
% Visualize results
figure('WindowState','maximized','Name','Displacement History')
hold on; grid on; box on
plot(T,U(:,1),'LineWidth',2)
xlabel('$T$ [s]','Interpreter','latex')
ylabel('$U$ [m]','Interpreter','latex')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex')
hold off