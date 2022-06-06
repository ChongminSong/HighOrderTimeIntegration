The Fortran code is called from Matlab and the data are exchanged as text files.

With MKL Pardiso, solving a complex equation (increasing M by 2) takes about twice the time of solving a real equation (increasing M by 1). The results are numerically identical with those of Matlab.

The Fortran code is written using the following
Visual Studio Community 2019
https://visualstudio.microsoft.com/vs/

Intel oneAPI Base toolkit
https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html?operatingsystem=window&distributions=webdownload

Intel oneAPI HPC toolkit
https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit/download.html

If you like, you can install them and I will pass you the code.

---------------------------------------------

There is also a Matlab wrapper file calling the Fortran code. You need to change the folder of the executable file and output text file to match your setup. You may test the executable without using the source code.

If you encounter an error message about a missing libXXX file, please add the following to the PATH of windows: 
C:\Program Files (x86)\Intel\oneAPI\compiler\2021.1.1\windows\redist\intel64_win\compiler

Regards,
Chongmin
