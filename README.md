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

## 4. Algorithm Review

Here we explain the main steps of the Dbar algorithm with either direct numerical implementation or with theoretical motivation.

### Dirichlet-to-Neumann map

The first step in EIT is to apply currents on a finite set of electrodes around the boundary and measure the corresponding voltages obtained thereafter. Let **L** be the number of electrodes, then we can apply at most **L-1** linearly independent current patterns (by current pattern we designate a vector of lenght **L** that describes the values of current applied in the electrodes). We designate by **C<sup>n</sup>** the n-th current pattern. There are a few choices that can be considered for current patterns, but for now we fix our focus in two. For **n=1,..., L-1** we have:
- Adjacent current patterns: 
<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%5Chspace%7B3cm%7D%20C%5En_l%20%3D%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20M%2C%5Cquad%20%5C%2C%20n%3Dl%20%5C%5C%20-M%2C%20%5C%3B%20n%3Dl&plus;1%20%5C%5C%20%5Cquad%5C%3B0%2C%20%5C%3B%20%5Ctext%7B%20otherwise%20%7D%20%5Cend%7Bmatrix%7D%5Cright." /> </p>

- Trignometric current patterns: 
<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%5C%2CC%5En_l%20%3D%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20M%5C%2C%5Ctext%7Bcos%7D%5Cleft%28n%5Ctheta_l%5Cright%29%2C%5Cquad%5Cquad%5Cquad%20%5C%2C%20n%3D1%2C...%2C%5Cfrac%7BL%7D%7B2%7D%20%5C%5C%20%5C%5C%20%5Cquad%20M%5C%2C%5Ctext%7Bsin%7D%5Cleft%28%5Cleft%28n-L/2%5Cright%29%5Ctheta_l%5Cright%29%2C%20%5Cquad%20n%3D%5Cfrac%7BL%7D%7B2%7D&plus;1%2C%20...%2C%20L-1%20%5Cend%7Bmatrix%7D%5Cright." /> </p>

The corresponding voltages are denoted by **V<sup>n</sup>** and each of them is a vector of dimension **L** (corresponding to the set of electrodes). One of the required assumptions is that for each current pattern the sum of the measured voltages equals zero.
We assume that an EIT device provides the data in this manner. 

For our algorithm we create a mapping from this data: Dirichlet-to-Neumman map, which establishes a relation between the Dirichlet boundary conditions (voltages) and the Neumann boundary conditions (currents) in conductivity model. 

**How to obtain the Dirichlet-to-Neumann map (DtoN map)?**

1. Normalize the measurements:

     **c<sup>n</sup>= C<sup>n</sup>/||C<sup>n</sup>||<sub>2</sub>** and **v<sup>n</sup>= V<sup>n</sup>/||C<sup>n</sup>||<sub>2</sub>**
     
2. Define the Neumann-to-Dirichlet map (inverse of DtoN map), which is a L-1 x L-1 matrix:

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%5C%2CR_%7B%5Cgamma%7D%28i%2C%20j%29%20%3D%20%28t%5Ei%2C%20v%5Ej%29_L%20%3D%20%5Csum_%7Bl%3D1%7D%5E%7BL%7D%20t%5Ei_lv%5Ej_l" /> </p>

3. Compute the matrix approximation of Dirichlet-to-Neumann map for electrodes of Area **A** and a body of radius **r**:

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5C%5C%20%5C%2CL_%7B%5Cgamma%7D%20%3D%20%5Cfrac%7BA%7D%7Br%7D%20%5Cleft%28R_%7B%5Cgamma%7D%20%5Cright%20%29%5E%7B-1%7D%20%5C%5C" /> </p>


This mapping is the matrix approximation of the continuum operator for a body with conductivity &gamma;.

Furthermore, the Dbar algorthm requires either the DtoN map of a body with the same geometry but with conductivity 1 or a reference DtoN map that corresponds to the same body with a different physiology (e.g. expiration and inspiration). This last case is designated by **tdEIT** and is very useful when monitoring health conditions that change with time like air ventilation in the lungs. With this in mind, the ***read_data class*** was constructed to compute one DtoN map from one set of current patterns and corresponding measured voltages. 

### Scattering Transform

 The second step of the algorithm is to compute the scattering transform. This transform is of non-physical nature, i.e. has no physical explanation, but is tightly connect with the the conductivity and can be determined by the Dirichlet-to-Neumann map. It is given in the new spectral parameter **k** in the following manner:
 
 <p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Ctextbf%7Bt%7D%28k%29%20%3D%20%5Cint_%7B%5Cpartial%5COmega%7D%20e%5E%7Bi%5Cbar%7Bk%7D%5Cbar%7Bz%7D%7D%5Cleft%28%5CLambda_%7B%5Cgamma%7D%20-%20%5CLambda_%7B1%7D%20%5Cright%20%29%5Cpsi%28z%2C%20k%29%5C%2C%20ds%28z%29%20%3D%20%5Cint_%7B%5Cmathbb%7BR%7D%5E2%7D%20e%5E%7Bi%5Cbar%7Bk%7D%5Cbar%7Bz%7D%7Dq%28z%29%5Cpsi%28z%2C%20k%29%20%5C%2Cdxdy%5C%2C" /> </p>

Now due to the ill-posedness of the inverse problem, the noise in the measurements highly affects the values of **t** for large **|k|**. The regularization strategy proven by Siltanen shows that we only need to compute the following approximation for the scattering transform:
<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Ctextbf%7Bt%7D%5ER%28k%29%20%3D%20%5Cleft%5C%7B%20%5Cbegin%7Bmatrix%7D%20%5Cint_%7B%5Cpartial%5COmega%7D%20e%5E%7Bi%5Cbar%7Bk%7D%5Cbar%7Bz%7D%7D%5Cleft%28%5CLambda_%7B%5Cgamma%7D%20-%20%5CLambda_%7B1%7D%20%5Cright%20%29%5Cpsi%28z%2C%20k%29%5C%2C%20ds%28z%29%2C%20%5Ctext%7B%20for%20%7D%20%7Ck%7C%3CR%20%5C%5C%20%5C%5C%200%2C%20%5Cquad%20%5Cquad%20%5Cquad%5Cquad%5Cquad%5Cquad%5Cquad%5Cquad%5Cquad%5Cquad%5Cquad%5Cquad%5Ctext%7B%20otherwise%7D%20%5Cend%7Bmatrix%7D%5Cright." /> </p>

The beauty of this regularization strategy lies in a finite discretization of **k** inside a disk of radius **R** being required. The choice of **R** dependends in an estimate of the noise present in the measurements.

To completely understand the computation of the scattering transform we need to obtain **&psi;** at the boundary. Nachman's was able to prove that the values of &psi; at the boundary can be determined by the solution of a boundary integral equation for **k &in; &#8450; \ {0}**:
<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cleft.%5Cpsi%28%5Ccdot%2C%20k%29%5Cright%7C_%7B%5Cpartial%5COmega%7D%20%3D%20%5Cleft.%20e%5E%7Bikz%7D%5Cright%7C_%7B%5Cpartial%5COmega%7D%20-%20S_%7Bk%7D%28%5CLambda_%7B%5Cgamma%7D-%7B%5CLambda_1%7D%29%5Cpsi%28%5Ccdot%2C%20k%29" /> </p>

However, this equation takes a lot of time to solve, thus for speedup purposes we choose an approximate solution. We can have two approaches:
- Solve the boundary integral equation with **S<sub>0</sub>** substituting **S<sub>k</sub>** and with the following definition:
<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cleft%28S_0%20f%20%5Cright%20%29%28x%29%20%3D%20%5Cint_%7B%5Cpartial%5COmega%7D%20G_0%28x%20-%20y%29f%28y%29%5C%2C%20d%5Csigma%28y%29%2C%20%5Ctext%7B%20where%20%7D%20G_0%28x%29%20%3D%20-%5Cfrac%7B1%7D%7B2%5Cpi%7D%5Ctext%7B%20log%7D%5C%2C%7Cx%7C." /> </p>

This approach is present in [2].

- Take **&psi;&sime; e<sup> ikz</sup>**, since this is the leading behavior and do not solve any integral equation.
This is the fastest way, but its triviality leads to more errors.

Now the last step is to solve the Dbar equation, which gives the name to the algorithm.

### D-bar equation

Now, having computed the scattering transform or an approximation of it, we are able to solve the following PDE in the **k** variable:
<p align="center">
 <img src="https://latex.codecogs.com/png.latex?%5Ctext%7BFor%20all%20%7D%20z%5Cin%5COmega%2C%20%5Ctext%7Bsolve%20the%20D-bar%20equation%3A%20%7D%5Cquad%5Cquad%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%5Cbar%7Bk%7D%7D%20%5Cmu%28z%2C%20k%29%20%3D%20%5Cfrac%7B%5Ctextbf%7Bt%7D%28k%29%7D%7B4%5Cpi%5Cbar%7Bk%7D%7De%5E%7B-2i%5Ctext%7BRe%7D%28kz%29%7D%5Coverline%7B%5Cmu%28z%2Ck%29%7D" /> </p>
 
Afterwards, we can determine the conductivity by:
<p align="center">
 <img src="https://latex.codecogs.com/png.latex?%5Csigma%28z%29%20%3D%20%5Cleft%28%5C%2C%5Cmu%28z%2C%200%29%5C%2C%5Cright%29%5E2"/> </p>.
 
 By substituting the scattering transform by its approximations and solving the equation with respect to it we obtain respective approximations to the conductivity.
 
 ## 5. Numerical Implementation
 
 Above, we present an immediate numerical implementation of the Dirichlet-to-Neumann map that trivially translates to code and a theoretical implementation of the regularization strategy for the D-bar and corresponding determination of the conductivity. Accordingly, we present the numerical implementation for both versions of the scattering transform and the solution to the Dbar equation.
 
 ### t-exp Scattering transform:
 
 The t-exp scattering transform corresponds to take the exponential approximation of &psi;, which has the following form:
 
 <p align="center">
  <img src="https://latex.codecogs.com/png.latex?%5Ctextbf%7Bt%7D%5E%7B%5Ctext%7Bexp%7D%7D%28k%29%20%3D%20%5Cint_%7B%5Cpartial%5COmega%7D%20e%5E%7Bi%5Cbar%7Bk%7D%5Cbar%7Bz%7D%7D%5Cleft%28%20%5CLambda_%7B%5Cgamma%7D%20-%20%5CLambda_1%20%5Cright%20%29e%5E%7Bikz%7D%5C%2C%20ds%28z%29"/> </p>
  
Let us consider that we have equally distributed electrodes around the boundary. By our assumption of &Omega; being a disc we define, without loss of generality, the positions of the electrodes to be <img src="https://render.githubusercontent.com/render/math?math=\vec{z}=(z^1,\, z^2, \, ...,\, z^L)" > where <img src="https://render.githubusercontent.com/render/math?math=z^j=(2\pi j)/L" >.

From this, we can take an approximation of the exponentials in these set points as an expansion in the orthonormal basis of current pattern (**L-1** degrees of freedom):

<p align="center">
 <img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20e%5E%7Bik%5Cvec%7Bz%7D%7D%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7BL-1%7D%20a_j%28k%29c%5Ej%20%5Cquad%20%5Cquad%20%5Cquad%20%5Ctext%7B%20and%20%7D%20%5Cquad%20%5Cquad%20%5Cquad%20e%5E%7Bi%5Cbar%7Bk%7D%5Cbar%7B%5Cvec%7Bz%7D%7D%7D%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7BL-1%7D%20b_j%28k%29c%5Ej" />. </p>
 
By immediate substitution and approximation of the integral at the electrode points this leads too:

<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20%5Ctextbf%7Bt%7D%5E%7B%5Ctext%7Bexp%7D%7D%28k%29%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7BL-1%7D%5Csum_%7Bi%3D1%7D%5E%7BL-1%7D%20a_i%28k%29b_j%28k%29%20%28c%5Ej%2C%20%5Cleft%28%5CLambda_%7B%5Cgamma%7D%20-%20%5CLambda_1%20%5Cright%20%29c%5Ei%29%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7BL-1%7D%5Csum_%7Bi%3D1%7D%5E%7BL-1%7D%20a_i%28k%29b_j%28k%29%20%5Cleft%28L_%7B%5Cgamma%7D%28j%2C%20i%29%20-%20L_1%28j%2C%20i%29%20%5Cright%29" />. </p>


### t-B Scattering Transform:

Here the first step is to solve the boundary integral equation with **S<sub>0</sub>**. For such we consider an initial vector approximation **&psi;**<sup>B</sup> of the solution of the (continuous) boundary integral equation and take the expansion in terms of the orthonormal current patterns:
<p align="center">
 <img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20%5Cpsi%5EB%28k%29%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7BL-1%7D%20p_j%28k%29%5C%2Cc%5Ej%20%3D%20%5Cleft%5Bc%5E1%5C%2C%20c%5E2%5C%2C%20...%20%5C%2Cc%5E%7BL-1%7D%20%5Cright%20%5D%5Ctextbf%7Bp%7D_k" />. </p> 
 
 Further, we also take a **LxL** matrix approximation of **G<sub>0</sub>** that defines the **S<sub>0</sub>** operator as:
 
 <p align="center">
  <img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20%5Ctextbf%7BG%7D_0%28i%2C%20j%29%20%3D%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20-%5Cfrac%7B1%7D%7B2%5Cpi%7D%5Ctext%7Blog%7D%20%5C%2C%7Cz%5Ei%20-%20z%5Ej%7C%2C%20%5Cquad%5Ctext%7B%20if%20%7D%20i%5Cneq%20j%20%5C%5C%20%5C%5C%200%5Cquad%5Cquad%5Cquad%5Cquad%5Cquad%5Cquad%5C%2C%2C%5Cquad%5C%2C%20%5Ctext%7B%20if%20%7D%20i%3Dj%20%5Cend%7Bmatrix%7D%5Cright." />. </p>
  
 Through this matrix approximation, we can translate the boundary integral equation into a linear system of equations. Defining **&Phi;** to be the current pattern matrix (**L x L-1**) leads us too:
 
 <p align="center">
  <img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20%5CPhi%5C%2C%5Ctextbf%7Bp%7D_k%20%3D%20e%5E%7Bik%5Cvec%7Bz%7D%7D%20-%20%5Ctextbf%7BG%7D_0%5Cleft%28L_%7B%5Cgamma%7D%20-%20L_1%20%5Cright%20%29%5CPhi%5C%2C%5Ctextbf%7Bp%7D_k%20%5C%5C%20%5C%5C%20%5CRightarrow%20%28%5CPhi%5ET%5CPhi%20&plus;%20%5CPhi%5ET%5Ctextbf%7BG%7D_0%5Cleft%28L_%7B%5Cgamma%7D%20-%20L_1%20%5Cright%20%29%5CPhi%29%5Ctextbf%7Bp%7D_k%20%3D%20%5CPhi%5ETe%5E%7Bik%5Cvec%7Bz%7D%7D" />, </p>
  
  
To compute the scattering transform we need to solve the last system of equations for each **k** in a disk of radius **R**. Afterwards, we compute it by:

<p align="center">
  <img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20%5Ctextbf%7Bt%7D%5EB%28k%29%20%3D%20A%20%5Csum_%7Bl%3D1%7D%5E%7BL%7D%20%5Csum_%7Bn%3D1%7D%5E%7BL-1%7D%20e%5E%7Bi%5Cbar%7Bk%7D%5Cbar%7Bz%7D%5El%7D%5Ctextbf%7Bp%7D_k%5En%28%5CPhi%5Cdelta%20L%20%29%28l%2C%20n%29" />, </p>
  where &delta;L is the difference between L<sub>&gamma;</sub> and L<sub>1</sub>.
 
 
 A more in depth look at the details can be found in [1] with respect to the solution of the boundary integral equation with **S<sub>k</sub>**.
 
 
 ### D-bar Equation:
 
 
 
 
 
 








 




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
 2. Add the t^exp approximation of the Scattering transform.

## 6. Bibliography:

[1] Mueller, J. L., & Siltanen, S. (2020). The d-bar method for electrical impedance tomography—demystified. Inverse problems, 36(9), 093001.
[2] El Arwadi, T., & Sayah, T. (2021). A new regularization of the D-bar method with complex conductivity. Complex Variables and Elliptic Equations, 66(5), 826-842.


