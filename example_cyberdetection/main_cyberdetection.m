%% Example of cyber-attack detection in power grids
% using functional observers and target state estimation in Ref. [1].
%
% References:
%
%   [1] A. N. Montanari, C. Duan, L. A. Aguirre, A. E. Motter. Functional
%       observability and target state estimation in large-scale networks.
%   [2] RD Zimmerman, CE Murillo-Sanchez, RJ Thomas, MATPOWER: Steady-State
%       Operations,993Planning, and Analysis Tools for Power Systems 
%       Research and Education.IEEE Transactions on Power Systems 26, 12?19 
%       (2011).
%   [3] T Nishikawa, AE Motter, Comparative analysis of existing models for
%       power-grid synchronization. New Journal of Physics 17, 015012 
%       (2015).

% Copyright (C) 2021  Arthur Montanari
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at
% your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% The full text of the GNU General Public License can be found in the 
% file license.txt.

%   Last modified by Arthur Montanari on 16/07/2021

clear all; close all; clc;
addpath('..');

%% Power grid model.
% The power grid dynamcis are modeled as a structure-preserving network of
% coupled Kuramoto oscillators (see Ref. [1], Eqs. 6-7 for details).

% Loads parameters (K,P,H,D,wref) used in Eqs. 6-7 (Ref. [1]). These
% parameters were computed based on the IEEE-118 power flow dataset
% (available in Ref. [2], MATPOWER case118.m) using the MATLAB toolbox 
% provided in Ref. [3] (using function SP_model.m).
load IEEE118pg.mat

% System dimensions
ng = 54;                        % number of generator buses
nl = 64;                        % number of load buses
N = 2*ng+nl;                    % number of Kuramoto oscillators
n = 3*ng+nl;                    % number of state variables

%% Continuous-time simulation of the power grid model in Eqs. 6-7.

% Time vector
t0 = 0.0;                       % initial time
tf = 20.0;                      % final time
dt = 0.0001;                    % integration step
tspan = t0:dt:tf;               % time span
T = length(tspan);              % number of data points

% Initial conditions for the generator phases and generator terminal phases
% are determined by the power flow solution (given by variable theta0 
% loaded from 1EEE118pg.mat). Generator frequencies and load buses phases 
% are initialized as zeros.
phi0 = [theta0; zeros(ng,1); theta0; zeros(nl,1)];

% Note that state vector phi is sorted as follows: [generators phases;
% generator frequencies; generator terminal phases; load phases].

% Power grid simulation using Runge-Kutta 4th order method.
[t,phi] = odeRK(@(t,theta)kuramotoPG_SP(t,theta,ng,nl,wref,K,H,P,D, ...
    gamma),[t0 dt tf], phi0');
xeq = phi(end,:)';              % equilibrium point/synchronization state

%% Defines the linear system matrix A 
% via linearization around synchronization state (operation point).
A = zeros(n,n);
A(1:ng,:) = [zeros(ng,ng) eye(ng) zeros(ng,ng+nl)];                         % matrix block [11 12 13]
A(ng+1:2*ng,ng+1:2*ng) = -D(1:ng)./(2*H).*eye(ng);                          % matrix block [   22   ]
thetaeq = [xeq(1:ng,1); xeq(2*ng+1:end,1)];
aux_ij = wref*K.*cos(ones(N,1)*thetaeq' - thetaeq*ones(1,N) + gamma).*(ones(N,N)./([2*H; D(ng+1:end)]));
aux_ii = - wref*eye(N).*sum(K.*cos(ones(N,1)*thetaeq' - thetaeq*ones(1,N) + gamma)).*(ones(N,N)./([2*H; D(ng+1:end)]));
A(ng+1:end,1:ng) = aux_ij(:,1:ng) + aux_ii(:,1:ng);                         % matrix block [21      ; 31      ]; 
A(ng+1:end,2*ng+1:end) = aux_ij(:,ng+1:end) + aux_ii(:,ng+1:end);           % matrix block [      23;       33]; 

%% PMU placement and cyber-attack location

% PMUs are randomly placed on the generator terminals and load buses. 
N_PMU = ceil(ng/3);            % number of PMUs
PMU = datasample(1:(ng+nl),N_PMU,'Replace',false);
S = PMU + 2*ng;                % set of sensor states (state variables
                               % which are directly measured by PMUs)
                               
% A single (randomly chosen) transmitted measurement by a PMU placed in 
% a generator terminal node is attacked. This is also the "set of target 
% nodes" T desired to be estimated by the functional observers.
TS = datasample(S(S<=3*ng),1,'Replace',false);
r = length(TS);
F = sparse(1:r, TS, 1, r, n);  % functional matrix (size rxn)

% Removes attacked nodes from the set of sensor nodes, i.e., attacked
% state variables are not used by a functional observer.
S = S(S~=TS);

% The attacked measurement is replaced by a false data copied from another
% neighboring measurement.
G = graph((K+K')/2);           % power grid (undirected) graph
Neigh = neighbors(G,TS - ng);  % finds neighbors of node TS-ng (which
                               % corresponds to state TS)
Neigh = Neigh(Neigh>ng) + ng;  % fake data comes only from loads/terminals
Fake = datasample(Neigh,1,'Replace',false);  % picks a random neighbor
                                             % to copy measurements for
                                             % false data injection

% Plots power grid topology, highlight generator nodes, generator
% terminals, load buses, sensor nodes (PMUs) and the cyber-attacked PMU.
figure(1); 
h = plot(G,'Layout','force','Iterations',40);
highlight(h,[1:1:ng],'NodeColor',0.1*[1 1 1]);    % plots generator nodes
highlight(h,[1:1:ng],'MarkerSize',6);
highlight(h,[ng+1:1:N],'NodeColor',0.55*[1 1 1]); % plots gen. terminals
highlight(h,[ng+1:1:2*ng],'MarkerSize',5);
highlight(h,[2*ng+1:1:N],'Marker','s');           % plots load nodes
highlight(h,[2*ng+1:1:N],'MarkerSize',6);
highlight(h,S-ng,'NodeColor',[0.011, 0.086, 0.827]);
highlight(h,S-ng,'MarkerSize',7);      % highlights PMUs
highlight(h,TS-ng,'NodeColor',[0.788, 0.149, 0.156]);
highlight(h,TS-ng,'MarkerSize',9);     % highlights cyber-attack
h.EdgeColor = [0.619 0.619 0.619];
h.LineWidth = 2;

%% System simulation during a cyber-attack as described in Ref. [1].

% Parameters
t0 = 0;                        % initial time
tf = 5;                        % final simulation time
tpert = 1;                     % time instant when system is perturbed
dt = 0.0001;                   % time step
tspan = t0:dt:tf;              % time span
T = length(tspan);             % number of data points

% System operates in steady-state until t=1s 
[~,x_ss] = odeRK(@(t,theta)kuramotoPG_SP(t,theta,ng,nl,wref,K,H,P,D, ...
    gamma),[t0 dt tpert-dt], xeq');

% System is hit by a small perturbation at time t = 1s.
x_pert = x_ss(end,:)' + [0.1*randn(ng,1); zeros(2*ng+nl,1)];
[~,x_af] = odeRK(@(t,theta)kuramotoPG_SP(t,theta,ng,nl,wref,K,H,P,D, ...
    gamma),[t0 dt tf-tpert], x_pert');

% System state evolution from t = 0s to t = 2s is given by
x = [x_ss; x_af]';

% Plots power grid dynamical evolution before and after perturbation
figure(2)
step = 100;

subplot(121)
plot(tspan(1:step:end),[x(1:ng,1:step:end); x(2*ng+1:end,1:step:end)])
title('Generators and loads phase evolution')
xlabel('t [s]'); ylabel('\phi_i [rad]')
xlim([t0 tf])

subplot(122)
plot(tspan(1:step:end),x(ng+1:2*ng,1:step:end))
title('Generators frequency evolution')
xlabel('t [s]'); ylabel('d\phi_i/dt [rad/s]')
xlim([t0 tf])

%% Target state estimation framework for cyber-attack detection.

Nobsv = 100;          % number of functional observers to be designed
                      % and implemented for the cyber-detection

% Dynamical system and observer simulation
for i = 1:Nobsv
    disp(['Iteration: ', num2str(i),'/',num2str(Nobsv),'.'])  
    
    % Output matrix: each functional observer is designed based on a random
    % set of half the total number of PMUs available.
    q = ceil(length(S)/2);
    Sdash = datasample(S,q,'Replace',false);
    C = sparse(1:q, Sdash, 1, q, n);
    
    % Minimum-order functional observer design
    F0 = find_F0(A,C,F);
    r0(i) = size(F0,1);
    [foN,foJ,foH,foD,foE] = functobsv_design(A,C,F0);
    
    % Inicialization (see Eq. 9, Ref. [1], for more details)
    y = zeros(q,T);                     % output vector             
    u = zeros(1,T);                     % input vector (null)
    zhat = zeros(r0(i),T);              % estimated target vector    
    w = zeros(size(foD,2),T);           % auxiliary vector
    
    % Target state estimation
    for k = 2:T
        % Output signal (measurements)
        y(:,k) = C*(x(:,k)-xeq) + 0.0*randn(q,1);    % no noise
        
        % Linear functional observer simulation (Eq. 9, Ref. [1])
        [~,waux] = odeRK(@(t,w)functobsv_sys(t,w,foN,foJ,foH, ...
            y(:,k-1),u(:,k-1)),[t0 dt dt],w(:,k-1)');
        w(:,k) = waux(end,:);
        zhat(:,k) = foD*w(:,k) + foE*y(:,k);
    end
    
    % Estimated target vector by the i-th functional observer
    zhatdash(i,:) = zhat(1:r,:) + F*xeq;  % shifts back coordinates to xeq
end

% Target states' true values
z = F*x;

%% Plot estimated target states versus true value and false-data injection.
figure(3);
step = 100;

hold on
plot(tspan(1:step:end),z(1:step:end),'-',...
    'Color',[0 0.4470 0.7410],'LineWidth',3.5)      % true value
plot(tspan(1:step:end),x(Fake,1:step:end) + (z(1) - x(Fake,1)),'-',...
    'Color',[0.788, 0.149, 0.156],'LineWidth',3)    % false data by attack
plot(tspan(1:step:end),zhatdash(1:1:end,1:step:end),'--',...
    'LineWidth',1)                                  % estimated values
xlabel('t [s]');
xlim([t0 tf])
