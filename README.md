The codes implement the high-order implicit time integration scheme found in the following two publications: 
1. Song, Ch., Eisenträger, S. and Zhang, X.R. “High-order implicit time integration scheme based on Padé expansions”, Computer Methods in Applied Mechanics and Engineering, Vol. 390, 114436. https://doi.org/10.1016/j.cma.2021.114436
2. Song, Ch., Zhang, X.R. and Eisenträger, S. “High-order implicit time integration scheme with controllable numerical dissipation based on mixed-order Padé expansions”, http://arxiv.org/abs/2206.04183 

The folder "PTIM_Matlab" contains a stand-alone implementation in Matlab. 

The folder "PTIM_Fortran" contains two subfolders "Fortran_Code" and "Matlab_Code". The Matlab functions prepare the inputs (stiffness matrix, mass matrix, force vectors, etc.) and output to temporary files. The Fortran program reads in these temporary files and performs the time integration. The folder names in the codes need to be modified to suit your own computer. 
