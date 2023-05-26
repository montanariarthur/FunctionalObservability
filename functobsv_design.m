function [N,J,H,D,E] = functobsv_design(A,C,F0,B,alpha,Qgain,Rgain)
% Designs the functional observer's matrices (N,J,H,D,E) in Eq. (9) in Ref.
% [1]. The design method, originally proposed in Ref. [2], is guaranteed 
% to provide a stable functional observer if the triple (A,C,F0) satisfies
% Darouach's conditions (4-5) in Ref. [1]. This code provides a MATLAB
% implementation of Algorithm 3.5.1 in Ref. [3].

% Inputs:
%   A          - dynamical matrix A (size nxn)
%   C          - output (measurement) matrix C (size qxn)
%   F          - functional matrix F0  (size r0xn)
%
% Optional:
%   B          - input (control) matrix; default = [0 0 ... 0]' (size 1xn)
%   alpha      - defines righ-most pole of matrix N; default = -100
%   Qgain      - Q matrix gain of a LQR for pole placement; default = 0.001
%   Rgain      - R matrix gain of a LQR for pole placement; default = 1
%
% Outputs:
%   N, J, H, D, E  - design matrices of a functional observer 
%                    (see Eq. 10, Ref. [1]).
%
% References:
%
%   [1] A. N. Montanari, C. Duan, L. A. Aguirre, A. E. Motter. Functional
%       observability and target state estimation in large-scale networks.
%       Proceedings of the National Academy of Sciences 119(1):e2113750119 (2022).
%   [2] M. Darouach. Existence and Design of Functional Observers for 
%       Linear Systems. IEEE Transactions on Automatic Control 45, 940?943 (2000).
%   [3] T. Hieu, F. Tyrone. Functional Observers for Dynamical Systems. 
%       Springer Berlin Heidelberg (2012).

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

%% Dimensions
n = size(A,1);
q = size(C,1);
r = size(F0,1);

% Input matrix
if nargin < 4
    B = sparse(ones(n,1));
end

%% Step 1: Compute sub-matrices A11, A12, A21, A22, F1, F2.
% Computes projection matrix
Cpinv = sparse(pinv(full(C)));  % Moore-Penrose pseudoinverse
Corth = sporth(spnull(C));      % sparse orthogonal basis of the null space
P = [Cpinv Corth];              % projection/transformation matrix
Pinv = inv(P);

% Matrix transformations
n_q = n - q;
Abar = Pinv*A*P; 
Fbar = F0*P;
A11 = Abar(1:q,1:q);
A12 = Abar(1:q,q+1:end);
A21 = Abar(q+1:end,1:q);
A22 = Abar(q+1:end,q+1:end);
F1 = Fbar(1:r,1:q);
F2 = Fbar(1:r,q+1:end);

%% Step 2: Compute N1 and N2 from Eq. 3.53b and 3.53c in Ref. [3].
F2ortho = sporth(spnull(F2));
Omega = A12*F2ortho;
Phi = -F2*A22*F2ortho;

% Eqs. 3.53b and 3.53c in Ref. [3]
Omegapinv = sparse(pinv(full(Omega)));
F2pinv = sparse(pinv(full(F2)));
N1 = (Phi*Omegapinv*A12 + F2*A22)*F2pinv;
N2 = (Omega*Omegapinv - speye(q))*A12*F2pinv;

%% Step 3: Determine Z such that N (Eq. 3.53a, Ref. [3]) is stable.
% Optimal pole placement method: LQR
if nargin < 5
    alpha = -100;
    Qgain = 1e-3;  
    Rgain = 1;
end
Q = Qgain*speye(r); 
R = Rgain*speye(size(N2',2));
Z = sparse(lqr(full(N1'+alpha*speye(r)),full(N2'),full(Q),full(R))');

N = N1 - Z*N2;

%% Step 4: Compute L1 from Eq. 3.52 in Ref. [3]. 
% Obtain L2 = F2, D = I_r.
% Obtain J and E from Eqs. 3.23 and 3.24 in Ref. [3].
L1 = Phi*Omegapinv + Z*(speye(q) - Omega*Omegapinv);
L2 = F2;
D = speye(r);
J = L1*A11 + L2*A21 - N*L1;
E = F1 - D*L1;

%% Step 5: Obtain H = LB, where L = [L1 L2]*inv(P).
L = [L1 L2]*Pinv;
H = L*B;

end
