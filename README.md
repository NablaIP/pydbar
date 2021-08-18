# pyDbar for Electrical Impedance Tomography

Thank you for the interest in pyDbar! This package implements the Dbar method for Electrical impedance tomography and follows up on the pyEIT presented in https://github.com/liubenyuan/pyEIT. 

## 1. Overview
The D-bar method can be used for Electrical Impedance Tomography in a 2D framework. This method is based on the theoretical reconstruction procedure of Nachman and a regularizing strategy that was established by Siltanen and provides a natural framework for numerical implementation There is an implementation in Matlab available at https://wiki.helsinki.fi/display/mathstatHenkilokunta/EIT+with+the+D-bar+method%3A+discontinuous+heart-and-lungs+phantom.

For now we will keep things simple and a proper documentation will follow.

## 2. Instalation

To install this package it is enough to use pip:
                      
                      pip install py-dbar.
 
 Moreover, dependencies are handled by pip.                    
                      
## 3. To start off:

Test cases are present in the tests folder. Soon, there will be a presentation here of the output.

## 4. Details on Implementation:

pyDbar is composed of three main classes: 

                      read_data, k_grid, dBar. 
                      
The class **read_data** and **k_grid** are independent of each other, while the **dBar** depends on both. To provide an overview of the functionality of all of them it is important to explain the theoretical reconstruction method.

### Dirichlet-to-Neumann map

The first step in EIT is to apply currents on a finite set of electrodes around the boundary and measure the corresponding voltages obtained thereafter. Let **L** be the number of electrodes, then we can apply at most **L-1** linearly independent current patterns (by current pattern we designate a vector of lenght **L** that describes the values of current applied in the electrodes). We designate by **T<sup>n</sup>** the n-th current pattern. There are a few choices that can be considered for current patterns, but we fix our focus for now in two. 
- Adjacent current patterns: T<sup> n, n</sup> = M and T<sup>n, n+1</sup> = -M and 0 otherwise, for n=1, ..., L-1
- Trignometric current patterns:
If T<sup> n, l</sup> = M cos(n$ \times (2 \pi l/L)$)

## 5. To ADD: 

 1. Support to simulate from pyEIT;
 2. Support of trignometric current patterns and pairwise current patterns;
 3. Support for GREIT and RL;
 4. Study U-Net to improve D-bar reconstruction.
