import numpy as np
import cmath
import math
# Contains Linear Operator and GMRES
import scipy.sparse.linalg as spla
from numpy import linalg as LG
import scipy.linalg as sla
from scipy.fft import fft2, ifft2
import pyamg
import matplotlib.pyplot as plt


pi = math.pi

class dBar:
    
    def __init__(self, R_z, m_z):
        self.Z = np.zeros( (2**m_z, 2**m_z), dtype=complex)
        self.load_mesh(1, m_z)
        
        self.sigma = np.zeros((2**m_z, 2**m_z))
    
    def load_mesh(self, R, m):
    
        N = int(pow(2, m))
        h = 2*R/N

        for l in range(N):
            for ll in range(N):
                    self.Z[l, ll] = complex(-R + l*h, -R + ll*h)
                
    
    
    def dBar(self, mu, k_grid, tK, zz):
    
        RHS = np.zeros((k_grid.N, k_grid.N), dtype=complex)

        N = len(k_grid.pos_x)

        for l in range(N):
                RHS[ k_grid.pos_x[l], k_grid.pos_y[l] ] = cmath.exp(-2j*( (k_grid.k[ k_grid.pos_x[l], k_grid.pos_y[l] ]*zz).real) )* tK.tK[k_grid.pos_x[l], k_grid.pos_y[l]]*complex(mu[l], -mu[l+N])


        F_RHS = fft2(RHS)

        for j in range(k_grid.N):
            for jj in range(k_grid.N):
                F_RHS[j, jj] = F_RHS[j, jj]*k_grid.FG[j, jj]

        RHS = ifft2(F_RHS)

        for l in range(N):
            mu[l] = mu[l] - (k_grid.h*k_grid.h)*RHS[ k_grid.pos_x[l], k_grid.pos_y[l]].real
            mu[l+N] = mu[l+N] - (k_grid.h*k_grid.h)*RHS[ k_grid.pos_x[l], k_grid.pos_y[l] ].imag


        return mu
    
  
    def solve(self, k_grid, tK):
        
        N = len(k_grid.pos_x)
    
        # Define the b and initial solution as the vectors of 1+0.j
        b = np.concatenate((np.ones(N), np.zeros(N)), axis=None)
        mu = np.concatenate((np.ones(N), np.zeros(N)), axis=None)

        
        zz = self.Z[0, 0]

        def Op(mu):
            return self.dBar(mu, k_grid, tK, zz)

        A = spla.LinearOperator((2*N,2*N), matvec=Op)

        n = self.Z.shape[0]
        
        for j in range(n):
            for jj in range(n):
                
                zz = self.Z[j, jj]

                mu, exitcode = pyamg.krylov.gmres(A, b, x0=mu, maxiter=5, orthog='mgs')
                #mu = EIT_GMRES(b, mu, 5, Kp, tK, Zz, FG)

                if(abs(zz) <= 1):
                    self.sigma[j, jj] = mu[k_grid.index]*mu[k_grid.index]-mu[k_grid.index+N]*mu[k_grid.index+N]


    def plot(self):
        
        Z_N = self.Z.shape[0]
        X = np.zeros(Z_N)
        Y = np.zeros(Z_N)
        h = 2/(Z_N)
        for i in range(Z_N):
            X[i] = - 1 + i*h
            Y[i] = - 1 + i*h
            
        sigma_x = np.zeros((Z_N, Z_N))

        for i in range(Z_N):
            for j in range(Z_N):
                sigma_x[j, i] = self.sigma[i, j]
            

        sigma_x = np.ma.masked_where(sigma_x==0, sigma_x)
        plt.pcolormesh(X, Y, sigma_x, cmap='RdBu')
        plt.colorbar()   
    
