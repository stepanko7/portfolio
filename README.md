# Portfolio


## Python implementaion of Zero Temperature String method in 2 and 3 dimensions

﻿In folder

    zero_temperature_string

one can find my implementation of the zero temperature string method,
which I have ported from Matlab.
The Matlab source code is presented at following link:

    https://cims.nyu.edu/~eve2/string.htm

During the porting process I have implemented bilinear and trilinear
interpolation methods:

    https://en.wikipedia.org/wiki/Bilinear_interpolation
    https://en.wikipedia.org/wiki/Trilinear_interpolation

The vectorized NumPy implementations are presented in

    bilinear_interpolate.py

and

    trilinear_interpolate.py

files.

Please, run

    python t0string.py

to see how the algorithm works in 2D and to examine the produced plots:

<img src="https://raw.githubusercontent.com/stepanko7/portfolio/master/zero_temperature_string/_final.png" width="500px" height="389px"/>
<img src="https://raw.githubusercontent.com/stepanko7/portfolio/master/zero_temperature_string/_fep1d.png" width="400px" height="311px"/>


You should have Python2.7, numpy, scipy and matplotlib installed.
If not, then run (in Linux):

    sudo pip install numpy
    sudo pip install scipy
    sudo pip install matplotlib

or install them manually in Windows.



## Visualization of methane adsorption isotherms and determination of optimal size

Folder

    adsorption_isotherms_visualization

contains an example Python script used to interpret and visualize
adsorption data (contained in "data" folder) obtained with
my implementation of the CDFT approach.
Please unpack data.7z or data.tar.gz
The script demonstrates utilization of main
NumPy methods such as stacking
of the 1D arrays to produce 2D array,
index selection and slicing.
Please run the script to produce the plots.
The data and figures were recently published in:

    S. Hlushak “Heat of adsorption, adsorption stress, and optimal storage of methane in slit and cylindrical carbon pores predicted by classical density functional theory” Phys. Chem. Chem. Phys. 20, 872--888 (2018)

<img src="https://raw.githubusercontent.com/stepanko7/portfolio/master/adsorption_isotherms_visualization/_deliver_dens_heat.png" width="500px" height="389px"/>
<img src="https://raw.githubusercontent.com/stepanko7/portfolio/master/adsorption_isotherms_visualization/_DRplot.png" width="500px" height="311px"/>



## Performing convolutions of 1D and 2D function in Fortran

Finally, in the folder

    convolutions_example

one could find a Fortran software which performs convolutions
of 1D and 2D tabulated "square" functions.
To compile and run the software in Linux, run:

    gfortran fftw_const.f90 convolutions.f90 ./libfftw3.a

and then

    ./a.out

The software will produce several files, which could be plotted
with the plot.py script:

    python plot.py

I have also supplied precompiled fftw library for linking in Linux.
You could also check the plot.py script for more Python code.

The convoluted functions and the result in the 1D-case
are presented in the following figure:

<img src="https://raw.githubusercontent.com/stepanko7/portfolio/master/convolutions_example/_convolution_1.png" width="320px" height="240px"/>



The convoluted functions in the 2D-case

<img src="https://raw.githubusercontent.com/stepanko7/portfolio/master/convolutions_example/_rho_2d.png" width="320px" height="240px"/>

<img src="https://raw.githubusercontent.com/stepanko7/portfolio/master/convolutions_example/_w3_2d.png" width="320px" height="240px"/>

and the corresponding result

<img src="https://raw.githubusercontent.com/stepanko7/portfolio/master/convolutions_example/_conv_2d.png" width="320px" height="240px"/>

are presented in the above three figures.
