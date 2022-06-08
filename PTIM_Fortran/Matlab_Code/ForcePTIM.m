function [ft] = ForcePTIM(order,signal,ndt,dt)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Chongmin Song & Sascha Eisentraeger (13/10/2021)
%
% SPDX-License-Identifier: MIT
%
%-------------------------------------------------------------------------%
% Input:
% order:            order of Pade expansion: = 11; 22; 33; 44; 55
% signal:           function handle for the excitation function
% ndt:          	number of time steps
% dt:               time steps size
%
% Output:
% ft:     Force matrix for the Pade time integration scheme; all values are
%         computed at t_i + dt/2 (middle of the time step); the matrix
%         contains the amplitude of the force and the derivatives
%         [amp damp/dt d^2amp/dt^2 d^(order/11)amp/dt^(order/11)];
%-------------------------------------------------------------------------%

switch order
    case 11
		% Interpolation points in the interval [-1,1]
        xi = gll_nodes(1);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
		% Calculate the force values within each time step
		t_vec1 = (s(1):1:ndt)'*dt;
        t_vec2 = (s(2):1:ndt)'*dt;
        amp_1 = signal(t_vec1);
		amp_2 = signal(t_vec2);
    	% Calculate the force vector of the time integration scheme:
        % The first column contains the value of the force at the middle of
        % the time step, while the second column holds the first derivative
        % of the force vector at the middle of the time step. --> The
        % forcing function is approximated by Lagrangian polynomials based
        % on GLL points; these functions are used to compute the values at
        % the midpoint.
        a1 = (s-0.5).^[0,1];
        A1 = (a1'*a1)\a1';
        ft = [amp_1(1:ndt),amp_2(1:ndt)]*A1';
    case 22
        % Interpolation points in the interval [-1,1]
        xi = gll_nodes(2);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
		t_vec1 = (s(1):1:ndt)'*dt;
        t_vec2 = (s(2):1:ndt)'*dt;
		t_vec3 = (s(3):1:ndt)'*dt;
        amp_1 = signal(t_vec1);
		amp_2 = signal(t_vec2);
		amp_3 = signal(t_vec3);
        % Calculate the force vector of the time integration scheme:
        % The first column contains the value of the force at the middle of
        % the time step, while the second column holds the first derivative
        % of the force vector at the middle of the time step, and the third
        % column is the second derivative of the force vector at the middle
        % of the time step. --> The forcing function is approximated by
        % Lagrangian polynomials based on GLL points; these functions are
        % used to compute the values at the midpoint.
        a2 = (s-0.5).^[0,1,2];
        A2 = (a2'*a2)\a2';
        ft = [amp_1(1:ndt),amp_2(1:ndt),amp_3(1:ndt)]*A2';
    case 33
        % Interpolation points in the interval [-1,1]
        xi = gll_nodes(3);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
		t_vec1 = (s(1):1:ndt)'*dt;
        t_vec2 = (s(2):1:ndt)'*dt;
		t_vec3 = (s(3):1:ndt)'*dt;
		t_vec4 = (s(4):1:ndt)'*dt;
        amp_1 = signal(t_vec1);
		amp_2 = signal(t_vec2);
		amp_3 = signal(t_vec3);
		amp_4 = signal(t_vec4);
        % Calculate the force vector of the time integration scheme:
        % The first column contains the value of the force at the middle of
        % the time step, while the second column holds the first derivative
        % of the force vector at the middle of the time step, and the third
        % column is the second derivative of the force vector at the middle
        % of the time step,... --> The forcing function is approximated by
        % Lagrangian polynomials based on GLL points; these functions are
        % used to compute the values at the midpoint.
        a3 = (s-0.5).^[0,1,2,3];
        A3 = (a3'*a3)\a3';
        ft = [amp_1(1:ndt),amp_2(1:ndt),amp_3(1:ndt),amp_4(1:ndt)]*A3';
    case 44
        % Interpolation points in the interval [-1,1]
        xi = gll_nodes(4);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
		t_vec1 = (s(1):1:ndt)'*dt;
        t_vec2 = (s(2):1:ndt)'*dt;
		t_vec3 = (s(3):1:ndt)'*dt;
		t_vec4 = (s(4):1:ndt)'*dt;
		t_vec5 = (s(5):1:ndt)'*dt;
        amp_1 = signal(t_vec1);
		amp_2 = signal(t_vec2);
		amp_3 = signal(t_vec3);
		amp_4 = signal(t_vec4);
		amp_5 = signal(t_vec5);
        % Calculate the force vector of the time integration scheme:
        % The first column contains the value of the force at the middle of
        % the time step, while the second column holds the first derivative
        % of the force vector at the middle of the time step, and the third
        % column is the second derivative of the force vector at the middle
        % of the time step,... --> The forcing function is approximated by
        % Lagrangian polynomials based on GLL points; these functions are
        % used to compute the values at the midpoint.
        a4 = (s-0.5).^[0,1,2,3,4];
        A4 = (a4'*a4)\a4';
        ft = [amp_1(1:ndt),amp_2(1:ndt),amp_3(1:ndt),amp_4(1:ndt),...
              amp_5(1:ndt)]*A4';
    case 55
        % Interpolation points in the interval [-1,1]
        xi = gll_nodes(5);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
		t_vec1 = (s(1):1:ndt)'*dt;
        t_vec2 = (s(2):1:ndt)'*dt;
		t_vec3 = (s(3):1:ndt)'*dt;
		t_vec4 = (s(4):1:ndt)'*dt;
		t_vec5 = (s(5):1:ndt)'*dt;
		t_vec6 = (s(6):1:ndt)'*dt;
        amp_1 = signal(t_vec1);
		amp_2 = signal(t_vec2);
		amp_3 = signal(t_vec3);
		amp_4 = signal(t_vec4);
		amp_5 = signal(t_vec5);
		amp_6 = signal(t_vec6);
        % Calculate the force vector of the time integration scheme:
        % The first column contains the value of the force at the middle of
        % the time step, while the second column holds the first derivative
        % of the force vector at the middle of the time step, and the third
        % column is the second derivative of the force vector at the middle
        % of the time step,... --> The forcing function is approximated by
        % Lagrangian polynomials based on GLL points; these functions are
        % used to compute the values at the midpoint.
        a5 = (s-0.5).^[0,1,2,3,4,5];
        A5 = (a5'*a5)\a5';
        ft = [amp_1(1:ndt),amp_2(1:ndt),amp_3(1:ndt),amp_4(1:ndt),...
              amp_5(1:ndt),amp_6(1:ndt)]*A5';
    otherwise
        % The variable order is incorrectly defined
        error("The value for order is not supported!");
end

end