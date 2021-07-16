function dx = SIRDnet(t,x,N,P,Adj,beta,gamma,eta)
% Represents the ODEs for a multi-group epidemiological model 
% (Eq. 9, Ref. [1]).
%
% Inputs:
%     t      -    Time
%     x      -    State vector (size 4*N)             
%     N      -    Number of populations/groups/cities/nodes
%     P      -    Population size (size N x 1)
%     Adj    -    Air transportation network (square matrix size NxN)
%     beta   -    Contact rate (size 2*ng+nl)
%     gamma  -    Inverse of infectious period (size N x 1)
%     eta    -    Case fatality rate (size N x 1)
%
% References:
%
%   [1] A. N. Montanari, C. Duan, L. A. Aguirre, A. E. Motter. Functional
%       observability and target state estimation in large-scale networks.

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

% State variables
S(:,1) = x(1:N)';                       % susceptible individuals S
I(:,1) = x(N+1:2*N)';                   % infected individuals I 
R(:,1) = x(2*N+1:3*N)';                 % recovered individuals R
D(:,1) = x(3*N+1:4*N)';                 % dead individuals D

% ODEs
dS = -beta.*S.*I./P - (S./P).*(Adj*ones(N,1)) + Adj'*(S./P);
dI = +beta.*S.*I./P - gamma.*I - (I./P).*(Adj*ones(N,1)) + Adj'*(I./P);
dR = (1-eta).*gamma.*I;
dD = eta.*gamma.*I;

% Returns
dx(1:N,1) = dS;
dx(N+1:2*N,1) = dI;
dx(2*N+1:3*N,1) = dR;
dx(3*N+1:4*N,1) = dD;

end