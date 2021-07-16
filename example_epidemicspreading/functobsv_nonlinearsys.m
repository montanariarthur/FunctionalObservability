function dw = functobsv_nonlinearsys(t,w,N,J,L,y,z,beta,n,P,F0)
% Represents the ODEs of a functional observer (Eq. S-24, Ref. [1]) 
% designed for the nonlinear multi-group epidemiological model (Eq. 9).
%
% Inputs:
%     t           -    Time
%     w           -    Auxiliary state vector             
%     N, J, H     -    Design matrices
%     y           -    Output vector
%     u           -    Input vector
%     beta        -    Contact rate
%     n           -    Number of groups/populations/nodes
%     P           -    Population size
%     F0          -    Functional matrix
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

% Nonlinear function f_1(v) = f_1(y,z), where v = [y; z] is an auxiliary
% vector, defined in S-25.
[Frow,Fcol] = find(F0);
for i = 1:n
    if sum(ismember(Fcol,i)) == 1 && sum(ismember(Fcol,i+n)) == 1
        aux1 = Frow(Fcol == i);
        aux2 = Frow(Fcol == i+n);
        v(i,1) = z(aux1)*z(aux2);
    else
        v(i,1) = 0;
    end
end

% Equation S-24.
dw = N*w'+J*y + L*[-beta.*v./P; +beta.*v./P; zeros(2*n,1)];

end