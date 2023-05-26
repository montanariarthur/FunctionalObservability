function dw = functobsv_sys(t,w,N,J,H,y,u)
% Represents the ODEs of a functional observer (see Eq. 10, Ref. [1] for
% details).
%
% Inputs:
%     t           -    Time
%     w           -    Auxiliary state vector             
%     N, J, H     -    Design matrices
%     y           -    Output vector
%     u           -    Input vector
%
% References:
%
%   [1] A. N. Montanari, C. Duan, L. A. Aguirre, A. E. Motter. Functional
%       observability and target state estimation in large-scale networks.
%       Proceedings of the National Academy of Sciences 119(1):e2113750119 (2022).

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

dw = N*w'+J*y+H*u;
end
