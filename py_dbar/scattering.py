import numpy as np
import cmath
import math

pi = math.pi

class scattering:
    
    def __init__(self, scat_type, k_grid, Now, Ref):
        self.tK = np.zeros((k_grid.N, k_grid.N), dtype=complex)
        
        self.load_scattering(scat_type, k_grid, Now, Ref)
        
        
    def load_scattering(self, scat_type, k_grid, Now, Ref):
        
        if scat_type=="partial":
            self.partial_scattering(Now, Ref, k_grid)
        elif scat_type=="exp":
            self.exp_scattering(Now,Ref, k_grid)
        else:
            print("Does not exist or we need to implement that version!")
            
            
        
    def exp_scattering(self, Now, Ref, k_grid):
        
        dt = (2*pi)/Now.L
        zt = np.exp(1j*np.arange(0, 2*pi, dt))
        ind0 = 1e-7
        
        for j in range(k_grid.N):
            for jj in range(k_grid.N):
                
                if( abs(k_grid.k[j, jj]) < k_grid.R and abs(k_grid.k[j, jj]) > ind0):
                
                    Ez = np.exp(1j*k_grid.k[j, jj]*zt)
                    conj_Ez = np.exp(1j * k_grid.k[j, jj].conjugate() * np.conjugate(zt))
                
                    ck, residuals_c, rank_c, s_c = np.linalg.lstsq(Now.Current, Ez, rcond=None)
                    dk, residuals_d, rank_d, s_d = np.linalg.lstsq(Now.Current, conj_Ez, rcond=None)
                
                    for l in range(Now.L-1):
                        for ll in range(Now.L-1):
                            self.tK[j, jj] = self.tK[j, jj] + (ck[l]*(Now.DNmap[ll, l]-Ref.DNmap[ll, l])*dk[ll])


                    self.tK[j, jj] = self.tK[j, jj]/(4*pi*(k_grid.k[j, jj].conjugate()))    
        
    
    
    def partial_scattering(self, Now, Ref, k_grid):
    
        G0 = np.zeros((Now.L, Now.L), dtype=complex)


        dt = (2*pi)/Now.L
        zt = np.exp(1j*np.arange(0, 2*pi, 2*pi/Now.L))

        for l in range(Now.L):
            for ll in range(Now.L):
                if l != ll:
                    G0[l, ll] = -(1/(2*pi))*cmath.log( abs( zt[l] - zt[ll] ) )

        

        dL = Now.DNmap - Ref.DNmap
        Phi = np.matmul(Now.Current.transpose(), Now.Current)
        PhidL = np.matmul(Now.Current, dL)

        ind0 = 1e-7

        M = Phi + np.matmul(Now.Current.transpose(),np.matmul(G0, PhidL))

        for j in range(k_grid.N):
            for jj in range(k_grid.N):
                if( abs(k_grid.k[j ,jj]) < k_grid.R and abs(k_grid.k[j, jj]) > ind0):

                    Ez = np.exp(1j*k_grid.k[j, jj]*zt)

                    psi_b, residuals, rank, s  = np.linalg.lstsq(M, np.matmul(Now.Current.transpose(),Ez), rcond=None)


                    for l in range(Now.L):
                        c = cmath.exp(1j*((k_grid.k[j, jj]*zt[l]).conjugate()))
                        for ll in range(Now.L-1):

                            self.tK[j, jj] = self.tK[j, jj] + c*PhidL[l, ll]*psi_b[ll]

                    self.tK[j, jj] = self.tK[j, jj]/(4*pi*(k_grid.k[j, jj].conjugate()))

    
            