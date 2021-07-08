import numpy as np
import cmath
import math
# Contains Linear Operator and GMRES
import scipy.sparse.linalg as spla
from numpy import linalg as LG
import scipy.linalg as sla
from scipy.fft import fft, ifft
import pyamg
import matplotlib.pyplot as plt


pi = math.pi

class dBar:
    
    def __init__(self, k_grid, Now, Ref, R_z, m_z):
        self.tK = np.zeros(k_grid.k.size, dtype=complex)
        self.Scat_B(Now, Ref, k_grid)
        self.Z = np.zeros(int(pow(2, 2*m_z)), dtype=complex)
        self.load_mesh(1, m_z)
        self.sigma = np.zeros(self.Z.size)
        
        
    def Scat_B(self, Now, Ref, k_grid):
    
        G0 = np.zeros((Now.L, Now.L), dtype=complex)


        dt = (2*pi)/Now.L
        zt = np.exp(1j*np.arange(0, 2*pi, 2*pi/Now.L))

        for l in range(Now.L):
            for ll in range(Now.L):
                if l != ll:
                    G0[l, ll] = -(1/(2*pi))*cmath.log( abs( zt[l] - zt[ll] ) )

        Ez = np.zeros(Now.L, dtype=complex)

        dL = Now.DNmap - Ref.DNmap
        Phi = np.matmul(Now.Current.transpose(),Now.Current)
        M = np.zeros((Now.L-1, Now.L-1), dtype=complex)
        PhidL = np.matmul(Now.Current,dL)

        ind0 = 1e-7

        M = Phi + np.matmul(Now.Current.transpose(),np.matmul(G0, PhidL))

        for l in range(k_grid.k.size):
                if( abs(k_grid.k[l]) < k_grid.R and abs(k_grid.k[l]) > ind0):

                    Ez = np.exp(1j*k_grid.k[l]*zt)

                    psi_b, residuals, rank, s  = np.linalg.lstsq(M, np.matmul(Now.Current.transpose(),Ez), rcond=None)


                    for n in range(Now.L):
                        c = cmath.exp(1j*((k_grid.k[l]*zt[n]).conjugate()))
                        for m in range(Now.L-1):

                            self.tK[l] = self.tK[l] + c*psi_b[m]*PhidL[n, m]

                    self.tK[l] = Now.AE*self.tK[l]/(4*pi*(k_grid.k[l].conjugate()))


    
    def dBar(self, mu, k_grid, zz):
    
        RHS = np.zeros((k_grid.N*k_grid.N), dtype=complex)

        N = len(k_grid.idx)

        for l in range(N):
            RHS[k_grid.idx[l]] = cmath.exp(-2j*( (k_grid.k[k_grid.idx[l]]*zz).real) )* self.tK[k_grid.idx[l]]*complex(mu[l], -mu[l+N])


        F_RHS = fft(RHS, norm="backward")

        F_RHS = np.multiply(F_RHS, k_grid.FG)

        RHS = ifft(F_RHS, norm="backward")


        for l in range(N):
            mu[l] = mu[l] - (k_grid.h*k_grid.h)*RHS[k_grid.idx[l]].real
            mu[l+N] = mu[l+N] - (k_grid.h*k_grid.h)*RHS[k_grid.idx[l]].imag


        return mu
    
    
    def load_mesh(self, R, m):
    
        N = int(pow(2, m))
        h = 2*R/N

        for l in range(N):
            if l%2==0:
                for ll in range(N):
                    self.Z[l*N+ ll] = complex(-R + l*h, -R + ll*h)
            else:
                for ll in range(N):
                    self.Z[l*N+ ll] = complex(-R + l*h, R - h*ll)
                
    
    
    # We are assuming that the data is in EIT_Data/Simulation/...
    def solve(self, k_grid):
        
        N = len(k_grid.idx)
    
        # Define the b and initial solution as the vectors of 1+0.j
        b = np.concatenate((np.ones(N), np.zeros(N)), axis=None)
        mu = np.concatenate((np.ones(N), np.zeros(N)), axis=None)

        
        zz = self.Z[0]

        def Op(mu):
            return self.dBar(mu, k_grid, zz)

        A = spla.LinearOperator((2*N,2*N), matvec=Op)


        for i in range(1, self.Z.size):
                zz = self.Z[i]

                mu, exitcode = pyamg.krylov.gmres(A, b, x0=mu, maxiter=10, orthog='householder')
                #mu = EIT_GMRES(b, mu, 5, Kp, tK, Zz, FG)

                if(abs(zz) <= 1):
                    self.sigma[i] = mu[k_grid.indx]*mu[k_grid.indx]-mu[k_grid.indx+N]*mu[k_grid.indx+N]


    def plot(self):
        
        Z_N = int(math.sqrt(self.Z.size))
        X = np.zeros(Z_N)
        Y = np.zeros(Z_N)
        h = 2/(Z_N)
        for i in range(Z_N):
            X[i] = - 1 + i*h
            Y[i] = - 1 + i*h
            
        sigma_x = np.zeros((Z_N, Z_N))

        for i in range(Z_N):
            if i%2==0:
                for j in range(Z_N):
                    sigma_x[i, j] = self.sigma[j+i*Z_N]
            else:
                for j in range(Z_N):
                    sigma_x[i, j] = self.sigma[(Z_N-1-j)+i*Z_N]

        sigma_x = np.ma.masked_where(sigma_x==0, sigma_x)
        plt.pcolormesh(X, Y, sigma_x, cmap='binary')
        plt.colorbar()   
    
