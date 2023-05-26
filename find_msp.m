function S = find_msp(A,T,Cand,k)  
% Finds the minimum set of sensor nodes required for the system functional
% observability with respect to a set of target nodes. This is a MATLAB
% implementation of Algorithm 1 in Ref. [1].
%
% Inputs:
%   A           - system matrix A of a dynamical system/network (nxn)
%   T           - array with target nodes' index
%   Cand        - array with candidate nodes' index for sensor placement
%   k           - maximum observable distance from sensors (default k = n)
%
% Outputs:
%   S           - array with sensor nodes' index
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

% Dimensions
n = size(A,1);                      % size of the system matrix
r = length(T);                      % cardinality of the target set
n_cands = length(Cand);             % number of candidate nodes
G = digraph((A~=0)');               % inference graph G(A)
    % Note that A is transposed in the code since MATLAB defines an 
    % adjacency matrix as A(i,j) = 1 if node i points to j, while we use 
    % the convention, in Ref. [1], that A(i,j) = 1 if node j points to i.
if nargin < 4
    k = n;
end
    
% Outer for-loop in Algorithm 1: Breadth-first search
R = sparse(n_cands,r);
for i = 1:r                         % for all x_i \in T
    % Find the set of reachable nodes R'_i in graph G(A) starting at node 
    % x_i \in T using a breadth-first search.        
    Ri_aux = bfsearch(G,T(i));
    
    % Inner for-loop in Algorithm 1. Matrix R defines the set of target
    % nodes that are reachable from the candidate nodes, i.e., R(i,j) = 1 
    % if target node x_i is reachable from node x_j, and 0 otherwise.
    for j = 1:n_cands
        if sum(Cand(j) == Ri_aux) == 1 && ...
                length(shortestpath(G,T(i),Cand(j))) <= k
            R(j,i) = 1;
        end
    end
end

% While-loop in Algorithm 1: Greedy algorithm
i = 0; 
cond_while = 0;
R0 = R;
while cond_while == 0
    i = i + 1;
    
    % Finds the number of reachable targets per candidate node (per column)
    R_sum = sum(R,2);
    % Finds the largest set of reachable targets among all candidate nodes
    R_max = max(R_sum);
    
    % The elements x_i with highest gain \Delta(x_i) are the elements with
    % the largest set of reachable targets:
    Delta_max = find(R_sum == R_max);
    S(i) = datasample(Delta_max,1);       % adds a single candidate node 
                                          % randomly to the sensors set 
                                          % in case there are multiple
                                          % candidates with a same gain
                                          
    % Checks while condition: is the union of all reachable targets by the
    % chosen set of sensor nodes S equal to the target set?
    if sum(ismember(T,T(find(sum(R0(S,:),1)>=1)))) == r
        cond_while = 1;
    elseif cond_while == 0 && i == n_cands
        cond_while = 1;
        disp(['WARNING: System cannot be functional observable for the provided set of candidate nodes.']);
    end
    
    % Removes the chosen element from the set of candidates for next
    % iterations
    R = R(:,~R(S(i),:) == 1);
end

% Output
S = Cand(S);

end
