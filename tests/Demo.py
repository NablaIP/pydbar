import py_dbar as pd
from py_dbar import read_data,k_grid, dBar 

# Read the boundary measurements and create the Dirichlet-to-Neumann map
# Frame to determine the conductivity: 
Now = read_data.read_data("Object1/", 1, 1, 16)
# Frame of Reference - Homogeneous
Ref = read_data.read_data("ObjectH/", 1, 1, 16)

# Establish a k_grid which corresponds to a spectral parameter
Kp = k_grid.k_grid(4, 6)

# Creates the model which corresponds to the Dbar equation in the spectral parameter at each z of a z_grid.
model = dBar.dBar(Kp, Now, Ref, 1., 6)

# Solves the problem for each z in the 2D plane, hence determining the conductivity at each z
model.solve(Kp)

# Plots the solution
model.plot()
