import numpy as np
import cmath
import math
from scipy.fft import fft, ifft


pi = math.pi

class k_grid:
    
    def __init__(self, R, m):
        self.R = R
        self.m = m
        self.s = 2.3*R
        self.h = 2*2.3*R/pow(2, m)
        self.N = pow(2, m)
        self.k = np.zeros(int(self.N*self.N), dtype=complex)
        self.FG = np.zeros(int(self.N*self.N), dtype=complex)
        self.k_gen()
        self.Green_FS()
        self.find()
            
    def k_gen(self):
            
            for l in range(self.N):
                for ll in range(self.N):
                    self.k[l*self.N + ll] = complex(-2.3*self.R + l*self.h, -2.3*self.R + ll*self.h)
                    
    def find(self):
    
        idx = []
        self.indx = -1

        for l in range(self.k.size):
            if abs(self.k[l]) < self.R:
                idx.append(l)
                if(abs(self.k[l]) < 1e-7):
                    self.indx = len(idx)-1

        self.idx = np.array(idx)
    
         
    
    def Green_FS(self):
        eps = self.s/10
        RR = (self.s - eps)/2
        ind0 = 1e-7

        FG = np.zeros(int(self.N*self.N), dtype=complex)

        for l in range(self.k.size):
            if(abs(self.k[l]) < self.s and abs(self.k[l])>ind0 ):

                FG[l] = 1./(pi*self.k[l])

                if(abs(self.k[l]) >= 2*RR):

                    FG[l] = FG[l]*(1. - (abs(self.k[l]) - 2*RR )/eps)


        G =  np.zeros(int(self.N*self.N), dtype=complex)
        G[:] = FG

        ss=0
        ll=0
        N = int(pow(2, self.m))


        for k in range(N):
            for j in range(N):
                ss = int(j + N/2) % N
                ll = int(k + N/2) % N
                FG[j+k*N] = G[ss+ll*N]


        self.FG = fft(FG, norm="backward")
