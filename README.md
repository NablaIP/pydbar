# pyDbar for Electrical Impedance Tomography

This package follows up on the pyEIT presented in https://github.com/liubenyuan/pyEIT. In here, I implement the D-bar method in the 2D dimensional case. This method is based on the theoretical reconstruction procedure of Nachman and a regularizing strategy, that also simplifies numerical implementantion, was initially provided by Siltanen. There is an implementation in Matlab available at https://wiki.helsinki.fi/display/mathstatHenkilokunta/EIT+with+the+D-bar+method%3A+discontinuous+heart-and-lungs+phantom.

For now we will keep things simple and a proper documentation will follow.

To test it out simply install through: pip install py-dbar.
A simple use case is provided in the tests folder together with Data.


## To ADD: 

 1. Support to simulate from pyEIT;
 2. Support of trignometric current patterns and pairwise current patterns;
 3. Support for GREIT and RL;
 4. Study U-Net to improve D-bar reconstruction.
