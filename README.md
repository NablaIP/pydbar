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
                      
The class **read_data** and **k_grid** are independent of each other, while the **dBar** depends on both. 

### Dirichlet-to-Neumann map

The first step in EIT is to apply currents on a finite set of electrodes around the boundary and measure the corresponding voltages obtained thereafter. Let **L** be the number of electrodes, then we can apply at most **L-1** linearly independent current patterns (by current pattern we designate a vector of lenght **L** that describes the values of current applied in the electrodes). We designate by **T<sup>n</sup>** the n-th current pattern. There are a few choices that can be considered for current patterns, but for now we fix our focus in two. For **n=1,..., L-1** we have:
- Adjacent current patterns: 
<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%5Chspace%7B3cm%7D%20T%5En_l%20%3D%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20M%2C%5Cquad%20%5C%2C%20n%3Dl%20%5C%5C%20-M%2C%20%5C%3B%20n%3Dl&plus;1%20%5C%5C%20%5Cquad%5C%3B0%2C%20%5C%3B%20%5Ctext%7B%20otherwise%20%7D%20%5Cend%7Bmatrix%7D%5Cright." /> </p>

- Trignometric current patterns: 
<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%5C%2CT%5En_l%20%3D%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20M%5C%2C%5Ctext%7Bcos%7D%5Cleft%28n%5Ctheta_l%5Cright%29%2C%5Cquad%5Cquad%5Cquad%20%5C%2C%20n%3D1%2C...%2C%5Cfrac%7BL%7D%7B2%7D%20%5C%5C%20%5C%5C%20%5Cquad%20M%5C%2C%5Ctext%7Bsin%7D%5Cleft%28%5Cleft%28n-L/2%5Cright%29%5Ctheta_l%5Cright%29%2C%20%5Cquad%20n%3D%5Cfrac%7BL%7D%7B2%7D&plus;1%2C%20...%2C%20L-1%20%5Cend%7Bmatrix%7D%5Cright." /> </p>

The corresponding voltages are denoted by **V<sup>n</sup>** and each of them is a vector of dimension **L** (corresponding to the set of electrodes). One of the required assumptions is that for each current pattern the sum of the measured voltages equals zero.
We assume that an EIT device provides the data in this manner. 

For our algorithm we create a mapping from this data: Dirichlet-to-Neumman map, which establishes a relation between the Dirichlet boundary conditions (voltages) and the Neumann boundary conditions (currents) in conductivity model. 

**How to obtain the Dirichlet-to-Neumann map (DtoN map)?**

1. Normalize the measurements:

     **t<sup>n</sup>= T<sup>n</sup>/||T<sup>n</sup>||<sub>2</sub>** and **v<sup>n</sup>= V<sup>n</sup>/||T<sup>n</sup>||<sub>2</sub>**
     
2. Define the Neumann-to-Dirichlet map (inverse of DtoN map), which is a L-1 x L-1 matrix:

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%5C%2CR_%7B%5Cgamma%7D%28i%2C%20j%29%20%3D%20%28t%5Ei%2C%20v%5Ej%29_L%20%3D%20%5Csum_%7Bl%3D1%7D%5E%7BL%7D%20t%5Ei_lv%5Ej_l" /> </p>

3. Compute the Dirichlet-to-Neumann map for electrodes of Area **A** and a body of radius **r**:

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5C%5C%20%5C%2C%5CLambda_%7B%5Cgamma%7D%20%3D%20%5Cfrac%7BA%7D%7Br%7D%20%5Cleft%28R_%7B%5Cgamma%7D%20%5Cright%20%29%5E%7B-1%7D%20%5C%5C" /> </p>


This mapping is the matrix approximation of the continuum operator for a body with conductivity &gamma;.

However, for the algorithm to run we either need the DtoN map of a body with the same geometry but with conductivity 1 or we need a reference DtoN map that corresponds to the same body with a different physiology (e.g. expiration and inspiration). This last case is designated by **tdEIT** and is very useful when monitoring health conditions that change with time like air ventilation in the lungs.

With this in mind, the ***read_data class*** was constructed to compute one DtoN map from one set of current patterns and corresponding measured voltages. 


**Input of read_data:** 
- ***str_Object***: name of folder that contains Current and Voltage files, for now we assume that all data is an EIT_Data folder.
- ***r***         : radius of the body (onlu holds for circular objects);
- ***AE***        : area of one electrode;
- ***L***         : number of electrodes;


**Output:**
- ***Current***: matrix of size L x L-1 of the normalized current patterns, each column represents one normalized current pattern;
- ***Voltage***: matrix of size L x L-1 of the normalized voltages, where the i-th column represents the voltage measured for the i-th current pattern;
- ***DNmap***  : matrix of size L-1 x L-1 that defines the Dirichlet-to-Neumann map for the measured data.








  


## 5. To ADD: 

 1. Change class read_data to just use r, AE and L as input since this will be the same for all frames of data when we do tdEIT -> This will decrease memory. Moreover, we probably do not need to also save the voltages in the class.

