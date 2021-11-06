from __future__ import division, absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt
import cmath
import math

import pyeit.mesh as mesh
from pyeit.mesh import quality
from pyeit.eit.utils import eit_scan_lines
from pyeit.eit.fem import Forward

class sim:

    def __init__(self, anomaly, L):
        self.anomaly = anomaly
        self.Current = np.zeros((L-1, L))
        self.Voltage = np.zeros((L-1, L))
        self.L = L


    def simulate(self):

        for j in range(self.L-1):
            self.Current[j, j] = -1
            self.Current[j, j+1] = 1
            
        """ 0. build mesh """
        mesh_obj, el_pos = mesh.create(self.L, h0=0.08)

        # extract node, element, alpha
        pts = mesh_obj['node']
        tri = mesh_obj['element']
        x, y = pts[:, 0], pts[:, 1]

        mesh_new = mesh.set_perm(mesh_obj, anomaly=self.anomaly, background=1.0)
        perm = mesh_new['perm']


        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        # draw mesh structure
        tpc=ax1.tripcolor(x, y, tri, np.real(perm),
        edgecolors='k', shading='flat', cmap="RdBu", alpha=0.5)
        fig.colorbar(tpc)
        


        """ 1. FEM forward simulations """
        # setup EIT scan conditions
        ex_mat = eit_scan_lines(self.L, 1)
        ex_mat = ex_mat[0:self.L-1]
        
        for i in range(self.L-1):
            ex_line = ex_mat[i]
            fwd = Forward(mesh_obj, el_pos)
            f, _= fwd.solve(ex_line, perm=perm)
            self.Voltage[i,:] = np.real(f[el_pos])

        
        
            
        
