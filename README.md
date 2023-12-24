# FortranThermoStats
This is a Fortran 90 program that allows one to read a molecule file (```input.f90```, in this case, is the water molecule) and calculate and print all the thermodynamic information of the molecule. In the ```parameters.f90``` there are the parameters of the main.

<p align="center">
 <img src="variae/H2O_Vibration.gif" alt="Water Vibrations">
</p>

This script is part of my duty to accomplish the Theoretical Chemistry Exams at University of Florence in 2017.

To see the theory behind it, in the ```variae``` folder there is the pdf ```thermo.pdf``` written by Joseph W. Ochterski which explains the Thermochemistry in Gaussian program, from which my program is inspired. In the above mentioned folder there are also the files of the water optimization and frequency calcolation in Gaussian09 (```water_opt_freq.gau``` is the input, ```water_opt_freq.log``` is the output).

### To run the script:
In my case, since I used an Apple M1 chip, I download ```gfortran``` compiler from the python ```anaconda``` package management in a dedicated environment.

To install the ```gfortran```:

1) ```conda create --name Fortran python=3.10.13```
1) ```conda activate Fortran```
1) ```conda install -c conda-forge gfortran```

To run the code you can use the MakeFile:

1) ```make compile``` to compile the module and main script
1) ```make run``` to run the code
1) ```make clean``` to clean all the executable and compile file

    Here below the comands run by the ```MakeFile```:

    1) Compile the module file

        ```gfortran -c parameters.f90```
    1) Compile the script 

        ```gfortran FortranThermoStats.f90 parameters.f90 -o FortranThermoStats.o```
    1) If needed, you have to allow the exectution of the script

        ```chmod +x FortranThermoStats.o```
    1) Run the code:

        ```./FortranThermoStats.o```
