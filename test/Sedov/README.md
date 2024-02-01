In this folder several test cases for the Sedov blast wave problem are present. The Sedov blastwave problem is a classic test problem for testing the 
ability of hydrodynamic solvers to handle strong shocks. In order to simulate a strong explosion, the pressure and the temperature at the cells in the origin
are set to very high values (in the order of 10e10) compared to the undisturbed region, which is set to ambient values. In order to simulate a point-like
explosion this high pressure region has to be limited to a few cells. A sensitivity to the number of cells is performed. When 1 cell is chosen as opposed
to 5 the simulation performs much better, especially when considering the pressure close to the origin, which is overestimated in the 5 cells case.



E = 85672.6 J
p = E/V * (gamma - 1) = E/(4/3 * pi * dr^3) = 8.18528×10¹² Pa
rho = 1.163 kg/m3 