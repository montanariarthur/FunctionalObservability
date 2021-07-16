function dphi = kuramotoPG_SP(t,phi,ng,nl,wref,K,H,P,D,gamma)
% Represents the ODEs for a structure-preserving model of a power grid as 
% a network of coupled Kuramoto oscillators (see Eqs. 6-7 in Ref. [1]).
%
% Inputs:
%     t     -    Time
%     phi   -    State vector (size 3*ng+nl)             
%     ng    -    Number of generator buses
%     nl    -    Number of load buses
%     K     -    Coupling matrix (square matrix size 2*ng+nl)
%     P     -    Power generation/consumption (size 2*ng+nl)
%     H     -    Inertia constants (size ng x 1)
%     D     -    Damping constants (size 2*ng+nl x 1)
%     gamma -    Phase shift constants (square matrix size 2*ng+nl)
%     wref  -    Reference frequency
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
N = 2*ng+nl;                           % number of Kuramoto oscillators
x1(:,1) = phi(1:ng);                   % generator phase variables
x2(:,1) = phi(ng+1:2*ng);              % generator frequency variables
x3(:,1) = phi(2*ng+1:2*ng+ng+nl);      % load phase variables

% Coupling term
aux = sum( sin( ones(N,1)*[x1; x3]' - [x1; x3]*ones(1,N) + gamma).*K ,2) ;

% Generators are modelled as 2nd-order Kuramoto oscillators:
dx1 = x2;                              % generators phases
dx2 = wref./(2*H(1:ng)).*( P(1:ng) - (D(1:ng)/wref).*x2 + aux(1:ng,1) );
                                       % generators frequencies
                                       
% Generators terminals and load buses are modelled as 1st-order Kuramoto
% oscillators:
dx3 = wref./D(ng+1:end).*( P(ng+1:end) + aux(ng+1:end,1) );
                                       % [generator terminals phases; 
                                       %  load buses phase]

% Returns
dphi = [dx1; dx2; dx3];

% Note that state vector phi is sorted as follows: [generators phases;
% generator frequencies; generator terminal phases; load phases].

end
