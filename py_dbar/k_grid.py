import numpy as np
import cmath
import math
from scipy.fft import fft2, ifft2, fftshift


pi = math.pi 

class k_grid:
    
    pos_x = []
    pos_y = []
    
    def __init__(self, R, m):
        self.R = R
        self.m = m
        self.s = 2.3*R
        self.h = 2*(2.3*R)/(2**m)
        self.N = 2**m
        self.index = -1
        self.k = np.zeros((self.N, self.N), dtype=complex)
        self.generate()
        self.FG = self.fund_sol()
        
    def generate(self):
        
        for j in range(self.N):
            for jj in range(self.N):
                
                self.k[j, jj] = complex(-self.s + j*self.h, -self.s + jj*self.h)
        
                if(abs(self.k[j, jj]) < self.R):
                    self.pos_x.append(j)
                    self.pos_y.append(jj)
                
                if(abs(self.k[j, jj]) < 1e-7):
                    self.index = len(self.pos_x)-1
                    
      
    
    def fund_sol(self):
        
        eps = self.s/10
        RR = (self.s-eps)/2
        i0 = 1e-7
        
        G = np.zeros((self.N, self.N), dtype=complex)
        
        for j in range(self.N):
            for jj in range(self.N):
                
                abs_k = abs(self.k[j, jj])
                
                if (abs_k < self.s and abs_k > i0):
                    G[j, jj] = 1/(self.k[j, jj]*pi)
                    
                    if(abs_k >= 2*RR):
                        G[j, jj] = G[j, jj]*(1 - (abs_k - 2*RR)/eps)
        
        return fft2(fftshift(G))
        