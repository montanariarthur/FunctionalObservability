# Functional Observability of Networks
Functional observability analysis, minimum sensor placement and minimum-order functional observer design on large-scale dynamical networks, including applications to power grids and epidemics. 

# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".


# Usage


The following codes are direct implementations of the minimum sensor placement algorithm and the minimum-order functional observer design proposed in Ref. [1].

- `find_msp` : Finds the minimum set of sensor nodes required for the system functional observability with respect to a set of target nodes. This is a MATLAB implementation of Algorithm 1 in Ref. [1].

- `find_F0` : Finds F0 with minimum-order such that Darouach's condition (4) in Ref. [1] is satisfied for a triple (A,C,F0). This is a MATLAB implementation of Algorithm 2 in Ref. [1].

- `functobsv_design` : Designs the functional observer's matrices (N,J,H,D,E) in Ref. [1, Equation (10)]. The design method is guaranteed to provide a stable functional observer if the triple (A,C,F0) satisfies Darouach's conditions (4-5) in Ref. [1].




The following examples illustrates numerical results of the algorithms described above in applications to large-scale complex networks, power grids and epidemics. For more details, see the full description in each file.

- `example_dynamicalnetwork` : Example of minimum sensor placement and minimum-order functional observer design for a random complex dynamical networks. This code examplifies how to:
    1. apply Algorithm 1 to determine the minimum set of sensor nodes for functional observability of a dynamical network with respect to a given set of target nodes; and
    2. apply Algorithm 2 to design a minimum-order functional observer.
Both algorithms are further described in Ref. [1] and codes `find_msp.m` and `find_F0.m`.

- `example_cyberdetection` : Example of cyber-attack detection in power grids using functional observers and target state estimation.

- `example_epidemicspreading` : Example of sensor placement and nonlinear functional observer design for prevalence rate of infection estimation in target populations during epidemics.




The above examples are dependent on the following codes:

- `odeRK` : Implements the numerical integration algorithm Runge Kutta 4th order.
- `functobsv_sys`,`functobsv_nonlinearsys` : Represents the ODEs of a linear functional observer (Ref. [1], Eq. (10)) and a nonlinear functional observer (Ref. [1], Eq. (S-24)), respectively.
- `kuramotoPG_SP`,`SIRDnet` : ODEs of the power grid model (Ref. [1], Eq. (6-7)) and the multi-group epidemiological model (Ref. [1], Eq. (9)), respectively.
- `spnull`,`sporth` :  Computes a sparse orthonormal basis for the null space of S and the range of S, respectively.



# References
1. Arthur N. Montanari, Chao Duan, Luis A. Aguirre, Adilson E. Motter (**2021**). Functional observability and target state estimation in large-scale networks.
