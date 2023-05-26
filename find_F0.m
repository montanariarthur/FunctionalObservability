function F0 = find_F0(A,C,F)
% Finds F_0 with minimum-order such that Darouach's condition (4) in
% Ref. [1] is satisfied for a triple (A,C,F_0). This is a MATLAB 
% implementation of Algorithm 2 in Ref. [1].
%
% Inputs:
%   A        - system matrix A of a dynamical system/network (size nxn)
%   C        - output (measurement) matrix C (size qxn), which defines
%              sensor nodes
%   F        - functional matrix F  (size rxn), which defines target nodes
%
% Outputs:
%   F0       - matrix F0 (size r0xn, where r0>=r), which defines the 
%              minimum order and structure of the functional observer
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
n = size(A,1);
q = size(C,1);
r = size(F,1);

% Initialization
F0 = F;
r0 = r;
[~,T] = find(F);     % defines set of target nodes
[~,S] = find(C);     % defines set of sensor nodes

% Removes columns of A corresponding to sensor and target nodes. 
% See Note 1 for further explanation below.
A(:,[T; S]) = 0;
G = C*A;             % G is only composed by rows of CA, since columns
                     % corresponding to nonzero elements in C,F are removed
F0A = F0*A;

% This is the while-loop in Algorithm 2 in Ref. [1]. This loop
% incrementally augments F0 until the structural (generic) rank of
% [C; C*A; F0] is equal to the structural rank of [C; C*A; F0; F0*A] --
% satisfying Darouach's condition (4).
cond_while = 0;
while cond_while == 0
    
    % Computes the structural rank of matrix G and [G; F0*A] using MATLAB's
    % function sprank. This is equivalent to building a corresponding
    % bipartite graph of a given matrix and finding the number of
    % right-matched nodes via the maximum matching algorithm.
    sprank_G = sprank(G);
    sprank_GF0A = sprank([G; F0A]);
    
    % If true, Darouach's condition (4) is satisfied for a triple (A,C,F0).
    if sprank_G == sprank_GF0A
        cond_while = 1;
    
    % If false, F0 needs to be augmented.
    else
        % Defines set M_2, where x_j\in M_2 if [FO*A]_ij is a nonzero entry
        [~,M2] = find(F0A);
        M2 = unique(M2);
        
        % Defines set of candidate nodes C = M2\M1
        Cand = ~ismember(M2,find(dmperm(G)));
        Cand = M2(find(Cand));
        if length(Cand) == 0
            [~,Cand] = find(dmperm(G));
            Cand = unique(Cand);
        end
        
        % Draws a random element of C
        T = [T; datasample(Cand,1,'Replace',false)];
        r0 = r0+1;
        
        % Update matrices. Column associated with the newly added target
        % node (drawn element from set C) is removed since their
        % corresponding node in the bipartite graph B is always
        % right-matched (see Note 1 below).
        A(:,T(end)) = 0;
        G(:,T(end)) = 0;
        F0A(:,T(end)) = 0;
        F0A(end+1,:) = sparse(1,T(end),1,1,n)*A;
    end
end

% Returns matrix F0
F0 = [F; sparse(1:r0-r,T(r+1:end),1,r0-r,n)];

end

%% Comments.

% Notation.
% Define G = [C; C*A; F0], and the bipartite graph B(V,X,E), where V is the
% set of nodes where each element corresponds to a row of G, V' is a set of
% nodes where each element also corresponds to a column of G, and
% (v_i,v'_j) is an undirected edge in E if G(i,j)~=0.

% Note 1.
% If node x_j is a sensor or target node, then there is some row in G
% (corresponding to C or F0, respectively) such that the j-th column entry
% is a nonzero entry. Since, by assumption, this is a unique nonzero entry 
% in this row, then the v'_j \in V' is a right-matched node, i.e.,
% v'_j \in M_1. Thus, sensor and target nodes are always right-matched
% nodes, i.e., it is not needed to check if they are right-matched or not.
% Therefore, to reduce the dimensionality of the problem, we remove columns
% of A (and, thereby, of C*A, F0*A and G) associated with sensor and target
% nodes.
