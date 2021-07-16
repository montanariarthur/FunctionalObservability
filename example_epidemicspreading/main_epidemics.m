%% Example of inference of prevalence rate of infection in epidemics
% using functional observers and target state estimation in Ref. [1].
%
% References:
%
%   [1] A. N. Montanari, C. Duan, L. A. Aguirre, A. E. Motter. Functional
%       observability and target state estimation in large-scale networks.
%   [2] TranStats database - T-100 Domestic Segment (All Carriers).
%       Accessed at transtats.bts.gov/DatabaseInfo.asp?DB_ID=111 on 
%       April, 2020.

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

%% Air transportation network of the domestic flights in the United States
% Loads data acquired from the TranStats database (Ref [2]). To reduce the
% dimensionality of the problem, only flights (edges) with more than 30
% people travelling from one city to another are considered. Cities (nodes)
% that do not have an influx and outflux of 30 people/day are excluded from
% the air transportation network. Thus, the US air transportation network
% is described by the connections between 81 cities.
load airnet.mat
%  Adj         -    Adjacency matrix describing the US air transp. network
%                   (note that A(i,j) > 0 if i->j, differently from the
%                    convention adopted in Ref. [1] where A(i,j)>0 if j->i)
%  N           -    Number of cities (nodes) considered in Adj
%  P           -    Population per city (node)
%  city_coord  -    City (node) geographical coordinates      

%% Epidemic spreading model

% Nodal dynamics are described by a SIRD model (Eq. 9, Ref. [1]) with the
% following parameters
N = size(Adj,2);                   % number of nodes/cities
n = 4*N;                           % number of states
eta = 0.01*ones(N,1);              % case fatality rate
gamma = 0.16*ones(N,1);            % inverse of infectious period
beta= 0.4 + 0.1*randn(N,1);        % contact rate

%% Structural model of the epidemic spreading model
% The structural model is given by the adjacency matrix of a graph
% (unweighted) where edges represent the linear and nonlinear interactions
% between states in the ODEs in Eq. 9, Ref. [1]. In other words, the
% unweighted adjacency matrix Astruct(i,j) ~= 0 if \dot x_i is a function
% of x_j. Note that the specific entries of Astruct(i,j) in the following
% code are irrelevant, only the structure of zero and nonzero entries is
% taken into account.

% Linear part
A = zeros(n,n);                                   % initialization
A(1:N,1:N) = - diag(Adj*ones(N,1)./P) + Adj'./P;  % S_i
A(N+1:2*N,N+1:2*N) = - diag(gamma) - diag(Adj*ones(N,1)./P) + Adj'./P;
                                                  % I_i
A(2*N+1:3*N,N+1:2*N) = diag((1-eta).*gamma);      % R_i
A(3*N+1:4*N,N+1:2*N) = diag(eta.*gamma);          % D_i

% Structural model (including nonlinear functions)
Astruct = A;                                      % includes linear part
Astruct(1:N,N+1:2*N) = -diag(beta);               % S_i
Astruct(1:N,1:N) = Astruct(1:N,1:N) - diag(beta); % S_i
Astruct(N+1:2*N,1:N) = +diag(beta);               % I_i
Astruct(N+1:2*N,N+1:2*N) = Astruct(N+1:2*N,N+1:2*N) + diag(beta);   % I_i
Astruct = sparse(Astruct);

%% Set of target cities is randomly selected in the US air network
nodesample = datasample(1:N,N,'Replace',false); % random sequence of cities
r = 15;                          % number of targets
targetcities = nodesample(1:r);  % random target cities
T = targetcities + N;            % target variables (i.e., target variables
                                 % corresponding to the infected states I_i
                                 % of target cities)
                                 
F = sparse(1:r, T, 1, r, n);     % functional matrix

%% Minimum sensor placement
Cand = find(~ismember(1:N,targetcities));  % targets cannot be sensor cands
Cand = Cand + 3*N;               % only dead states D_i can be candidates
S = find_msp(Astruct,T,Cand,3);  % minimum sensor placement
sensorcities = S - 3*N;          % sensor city node index
q = length(S);                   % number of sensor cities
C = sparse(1:q, S, 1, q, n);     % output matrix

% Checks if the system is structurally observable (Eq. 2)
if sprank(obsv(Astruct,C)) == n
    disp(['Triple (A,C,F) is structurally observable.'])
else
    disp(['Triple (A,C,F) is NOT structurally observable.'])
end

% Checks if the system is structurally functional observable (Eq. 3)
if sprank([obsv(Astruct,C); F]) == sprank(obsv(Astruct,C))
    disp(['Triple (A,C,F) is structurally functional observable.'])
else
    disp(['Triple (A,C,F) is NOT structurally functional observable.'])
end

% Plots network structure, highlight population/city nodes, target nodes 
% and sensor nodes.
figure(1);
Graph = graph((Adj+Adj')/2);
h = plot(Graph,'NodeLabel',[],'EdgeColor','black','EdgeAlpha',0.10,...
    'XData',city_coord(:,2),'YData',city_coord(:,1),...
    'LineWidth',5*Graph.Edges.Weight/max(Graph.Edges.Weight) );
hold on;
highlight(h,1:N,'NodeColor','black')
h.MarkerSize = 5*(log(P) - log(min(P))+0.01);
highlight(h,sensorcities,'NodeColor','blue')
highlight(h,targetcities,'NodeColor','red')
highlight(h,targetcities(ismember(targetcities,sensorcities)),'NodeColor',[0.3010 0.7450 0.9330])
hold off;

%% Functional observer design
% Determines the minimum-order functional observer structure
F0 = find_F0(sparse(Astruct),C,F);
r0 = size(F0,1);
[~,F0col] = find(F0);

% Functional observer design for a class of nonlinear systems in Eq. S-19.
% The following steps provides the design of matrices (N,J,L,D,E) for a
% functional observer, stated in Eq. (S-20), designed *especifically* for
% the nonlinear multi-group epidemiological model studied in Eq. 8. See
% Ref. [1] and Supplementary Materials for more details.

% Step 1: Nonlinear function decomposition
[Crow,Ccol] = find(C);
[Frow,Fcol] = find(F0);
W = zeros(n,n-q-r0);         
count = 1;
for i = 1:N
    if sum(Ccol == i) ~= 1 & sum(Fcol == i) ~= 1
        W(i,count) = -beta(i)/P(i);
        W(i+nb,count+1) = +beta(i)/P(i);
        count = count + 2;
    end
end

% Step 2: Lipschitz constant
lipschitz = max(beta./N);

% Step 3: Matrix transformations
A = full(A);    
C = full(C);    
F = full(F);    
F0 = full(F0);            
Cpinv = pinv(C);                % Moore-Penrose pseudoinverse
Corth = orth(null(C));          % orthogonal basis of the null space
foP = [Cpinv Corth];            % projection/transformation matrix
Pinv = inv(foP);

n_q = n - q;                
Abar = foP*A*Pinv;
Fbar = F0*foP;                
Wbar = foP*W;
A11 = Abar(1:q,1:q);        A12 = Abar(1:q,q+1:end);
A21 = Abar(q+1:end,1:q);    A22 = Abar(q+1:end,q+1:end);
F1 = Fbar(1:r0,1:q);        F2 = Fbar(1:r0,q+1:end);
W1 = Wbar(1:q,:);           W2 = Wbar(q+1:end,:);
P1 = Pinv(1:q,:);           P2 = Pinv(q+1:end,:);

F2orth = orth(null(F2));    F2pinv = pinv(F2);
Omega = [A12*F2orth W1];    Omegapinv = pinv(Omega);
Phi = -[F2*A22*F2orth F2*W2];
N1 = (Phi*Omegapinv*A12 + F2*A22)*F2pinv;
N2 = (Omega*Omegapinv - eye(q))*A12*F2pinv;
L1bar = Phi*Omegapinv*P1 + F2*P2;
L2bar = (Omega*Omegapinv - eye(q))*P1;

% % Step 4: Solves LMI in Eq. S-24

% Defines variables
setlmis([])
Q = lmivar(1,[size(N1,1) 1]);
G = lmivar(2,[size(N1,1) size(C,1)]);
beta1 = 1;                  beta2 = beta1;

% Defines problem
lmiterm([1 1 1 Q],1,N1);        lmiterm([1 1 1 Q],N1',1);
lmiterm([1 1 1 G],-1,N2);       lmiterm([1 1 1 -G],N2',-1);
lmiterm([1 1 1 0],(lipschitz^2)*(beta1+beta2)*eye(size(N1,1)));
lmiterm([1 1 2 Q],1,L1bar);     lmiterm([1 1 3 G],1,L2bar);
lmiterm([1 2 1 Q],L1bar',1);    lmiterm([1 2 2 0],-beta1*eye(n));
lmiterm([1 3 1 -G],L2bar',1);   lmiterm([1 3 3 0],-beta2*eye(n));

% Solves LMI
LMISYS = getlmis;
[tmin,xfeas] = feasp(LMISYS);
Q = dec2mat(LMISYS,xfeas,Q);
G = dec2mat(LMISYS,xfeas,G);

% Step 5: Computes auxiliary matrices
Z = inv(Q(r0,r0))*G;
T1 = Phi*Omegapinv + Z*(speye(q) - Omega*Omegapinv);
T2 = F2;

% Step 6: Designs functional observer matrices (N,J,L,D,E)
foN = N1 - Z*N2;
foD = eye(r0);
foL = L1bar - Z*L2bar;
foJ = T1*A11 + T2*A21 - foN*T1;
foE = F1 - foD*T1;

%% Continuous time-simulation of the epidemic spreading
% Time vector
t0 = 0.0;                       % initial time
tf = 200;                      % final time
dt = 0.1;                      % integration step
tspan = t0:dt:tf;               % time span

% Epidemic outbreak
cityoutbreak = 76;              % node (Miami) where the outbreak starts
I0 = 1000;                      % initial number of infected individuals
x0 = zeros(n,1);                % x(0) = [S(0); I(0); R(0); D(0)]
x0(1:N,1) = P;                  % S_i(0)
x0(cityoutbreak,1) = P(cityoutbreak) - I0; % S_i(0) at outbreak city
x0(N+cityoutbreak,1) = I0;      % I_i(0) at outbreak city

% Epidemic spreading simulation
[t,x] = odeRK(@(t,s)SIRDnet(t,s,N,P,Adj,beta,gamma,eta),[t0 dt tf], x0');
x = x';
z = F*x;              % true values of I_i per target city i
for j = 1:r
    [~,arg] = max(z(j,:));
    tp(1,j) = t(arg); % time instant of the epidemic peak per target city
end

% Prediction/estimation of the "true" time instant of the epidemic peak. 
% We simulate the epidemiological model in Eq. 8 with arbitrary initial 
% conditions (different from the above) to contrast results where the 
% prediction is provided by a free simulation run of the model (without a 
% functional observer) and where the prediction/estimation is provided by 
% a functional observer. Results are computed for the average of Nmc Monte
% Carlo runs.
Nmc = 1;                       % number of Monte Carlo runs
for i = 1:Nmc
    % Arbitrary initial conditions
    cityoutbreak_false = datasample(1:N,1); % a city is chosen randomly as
                                            % the outbreak starting point
    I0_false = 1;
    x0_false = zeros(n,1);
    x0_false(1:N,1) = P;                    % S_i(0) - false prediction
    x0_false(cityoutbreak_false,1) = P(cityoutbreak_false) - I0_false;
    x0_false(N+cityoutbreak_false,1) = I0_false;
    
    % Prediction of the epidemic spreading by a free simulation run of the
    % model with arbitrary conditions (without functional observer)
    [~,x_free] = odeRK(@(t,s)SIRDnet(t,s,N,P,Adj,beta,gamma,eta),...
        [t0 dt tf], x0_false');
    x_free = x_free';
    
    % Estimation of the epidemic spreading provided by a functional 
    % observer fed with measurements of the sensor cities
    y = C*x;                                    % output vector
    w = zeros(size(foD,2),length(tspan));       % initialization
    w(:,1) = inv(foD)*(F0*x0_false - foE*C*x0); % init. with false assump.
    zhat = zeros(r0,length(tspan));             % initialization
    for k = 2:length(tspan)
        [~,waux] = odeRK(@(t,w)functobsv_nonlinearsys(t,w,foN,foJ,foL,...
            y(:,k-1),zhat(:,k-1),beta,N,P,F0),[t0 dt dt],w(:,k-1)');
        w(:,k) = waux(end,:);                   % auxiliary vector
        zhat(:,k) = foD*w(:,k) + foE*y(:,k);    % estimated target states
    end
    
    % Results: Compares the time instant where the epidemic peak actually
    % happens in each city with the epidemic peaks predicted with and
    % without the use of a functional observer.    
    z_free = F*x_free;   % prediction with free simulation
    for j = 1:r
        % Free simulation prediction
        [~,arg] = max(z_free(j,:));
        tp_free(i,j) = t(arg);
        % Functional observer estimate
        [~,arg] = max(zhat(j,:));
        tp_fobsv(i,j) = t(arg);
    end
end

%% Plot results
disp(['Time instant of the epidemic peak per target city:'])
disp(['True value --- Free simulation prediction --- Estimation with functional observer'])
results_abs = [tp' mean(tp_free,1)' mean(tp_fobsv,1)']

disp(['Average error between the true time instant of the epidemic peak and the predictions per target city:'])
disp(['Free simulation prediction --- Estimation with functional observer'])
results_dif = [abs(results_abs(:,1) - results_abs(:,2)) abs(results_abs(:,1) - results_abs(:,3))]

disp(['Average error of the prediction methods over all target cities'])
disp(['Free simulation prediction --- Estimation with functional observer'])
results_avg = mean(results_dif,1)


% Shows the state evolution and corresponding predictions, only for the
% last Monte Carlo run
figure(2)
subplot(311); plot(t,z)
xlabel('t [days]'); xlim([t0 tf]);
title('Epidemic spreading')
subplot(312); plot(t,z_free)
xlabel('t [days]'); xlim([t0 tf]);
title('Prediction based on free simulation')
subplot(313); plot(t,zhat(1:r,:))
xlabel('t [days]'); xlim([t0 tf]);
title('Prediction with functional observer')
