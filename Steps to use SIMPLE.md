**Steps to Run SIMPLE and plot the data**

First of all you need to down load the code from the website <https://github.com/itpplasma/SIMPLE>

SIMPLE computes statistical losses of guiding-center orbits for particles of given mass, charge and energy from the volume of 3D magnetic configurations. Orbits are traced via a symplectic integrator that guarantees conservation of invariants of motion within fixed bounds over long integration periods. A classifier based on Poincarè plots allows for accelerated prediction for confinement of regular (non-chaotic) orbits. There are 3 kinds of classifiers used in the code, namely Mikowski dimensions, topological orbit classifier and parallel invariant classification. The details of these can be found in respective articles.

The main focus of SIMPLE is to make computation of fusion alpha particle losses fast enough to be used directly in stellarator optimization codes. For that reason currently 3D configurations with nested magnetic flux surfaces are supported based on a VMEC equilibrium file in NetCDF format, and orbits are computed without taking collisions into account.
## **Building**
The build system for SIMPLE is CMake.

Required libraries:

- NetCDF
- LAPACK/BLAS

Supported compilers:

- GNU Fortan
- Intel Fortran

Then run it by the following commands

cd /path/to/SIMPLE

mkdir build

cd build

cmake ..

make

This will produce libsimple.so and simple.x required to run the code. Now you should run the main executable simple.x but you must have the input file simple.in as well as the VMEC equilibrium file in the same folder (the file name has to be changed to wout.nc).

For example in my case:     

khanm@faepop28 build $ ./simple.x         (don’t miss the ./ in the front)

Notes:

So you can take a VMEC equilibrium of a stellarator, e.g., from the website wout\_LandremanPaul2021\_QA\_<https://github.com/>

However, Obviously, the file name has to be changed to wout.nc,   or change the file name

in the input file to whatever equilibrium file you want to analyze.

You will get a start.dat with the phase-space coordinates of the particles and a times\_lost.dat and confined\_fraction.dat respectively describing the times when the particles got lost, and how the confined fraction of trapped and passing particles evolves over time.

The main output is confined\_fraction.dat, containing four columns:

1. Physical time
1. Confined fraction of passing particles
1. Confined fraction of trapped particles
1. Total number of particles

The sum of 2. and 3. yields the overall confined fraction at each time.

In addition start.dat is either an input for given or an output for randomly generated initial conditions. Diagnostics for slow convergence of Newton iterations are written in fort.6601.

**.dat files in SIMPLE**

Input file simple.in

\* tcut: Cut time for Minkowski classification

\* fast\_class: Use new classifiers (2023), don't need tcut

\* class\_plot: exit after Minkowski classification to generate

"Brazilian flag" plot (Fig. 8 of Accelerated Methods paper)

Data columns in files:

start.dat

\1) zstart(1,ipart)=r (normalized toroidal flux s)

\2) zstart(2,ipart)=theta\_vmec

\3) zstart(3,ipart)=varphi\_vmec

\4) normalized velocity module z(4) = v / v\_0:

zstart(4,ipart)=1.d0

\5) starting pitch z(5)=v\_\parallel / v:

times\_lost.dat

\1) Particle index\. Corresponds to line number in start\.dat

\2) Time t\_loss [s] when the particle is lost\.

\* If never lost or classified as regular, maximum tracing time trace\_time is written.

\* If ignored due to contr\_pp, which defines deep passing region as confined, -1 is written.

\3) Trapping parameter trap\_par that is 1 for deeply trapped, 0 for tp boundary and negative for passing\. Eq\. (3\.1) in Accelerated Methods paper\.

Whenever trap\_par < contr\_pp, particle is not traced and counted as confined.

confined\_fraction.dat

\1) Time in s, according to number of recording timesteps ntimstep\.

\2) confpart\_pass: Number of confined passing particles / total number of particles

\3) confpart\_trap: Number of confined trapped particles / total number of particles

\4) total number of particles

Col 2 + Col 3 gives Number of confined particles / total number of particles

In the first line, the total fractions of trapped and passing particles are written (all confined in 1st step).


**Generating the Figures of**

[Journal of Plasma Physics,](https://www.cambridge.org/core/journals/journal-of-plasma-physics) [Volume 89,](https://www.cambridge.org/core/journals/journal-of-plasma-physics/volume/5F7DA81D13DF167B4FDC257A92A38489) [Issue 3,](https://www.cambridge.org/core/journals/journal-of-plasma-physics/issue/01A78E61D8EDB04F1FA38F9C274BB8E6) June (2023), 955890301

First, please go to the folder ISHW\_2022/data/CLASSIFICATION/  and open the file README.TXT, all the important information is available there. We just copy that here for your information

The 2D classification plots (in the JPP paper) are produced in sub-directories RUN\_CLASS contained in each directory of a given magnetic configuration. For this:

\1) Run simple\.x in RUN\_CLASS subdirectory;

(e.g. ISHW\_2022/data/CLASSIFICATION/QH\_Drevlak/RUN\_CLASS/)

\2) Run a post-processor, \./dens\_cl\_glob\.x (link to executable is contained in each subdirectory);

\3) Start Matlab and run there a script gl\.m \.

This results in two \*.png files and two \*.pdf files (the same as in \*.png) in RUN\_CLASS subdirectory with classification by J\_parallel and by ideal/non-ideal.

WARNINGS:

\1) Note that reference magnetic field (5 T in all cases) is enforced on a starting surface\. For Subbotin configurations, starting surface was at s=0\.6 where B\_00=5\.5 T in case B\_00=5 T on axis (computations were at 10% lower field)\. We keep the same s=0\.6 for all other configurations for consistency\.

\2) Configuration of Subbotin has been scaled to 1000 cubic meters of the plasma volume\. Configurations of Drevlak were taken as is so that plasma volume is 1900 cubic meters for QI, 1863 cubic meters for QH and 1900 cubic meters for QA (about 1\.9 times larger than for Subbotin configurations)\.

\3) All computations were using a new method for healing the axis\. For the QH configuration, this method makes a visible change: mode (m,n)=(1,20) is "healed" up to point 51 of 64 because of irregularity at the edge (forther in the core it is smooth)\. For this reason, results with old method for healing axis have been added to QH and QI of Drevlak, in sub-directories OLD contained in RUN\_CLASS\.

- Most of the plotting  (.py) files are located in the folder /ISHW\_2022/data/

**Figures 7-12 (a or b):**

1. Use the file   plot\_classification.py
1. See lines 49-51 for data files location, the data files are prompt1.dat, regular1.dat

and stochastic1.dat

the mentioned data files are generated as discussed above in green color.

1. Upon running the code in Visual Studio (interactive mode), we get the following plot (you can also run the code in Matlab as discussed in above green color text)




**Figures 7-12 (c):**

1. Use the file   plot\_loss\_over\_time.py
1. See line 36 for data file location, the data file is   energyslow\_aver.dat

the mentioned data file is generated by a F90 code   energy\_loss\_fraction.f90

1. Upon running the code in Visual Studio (interactive mode), we get the following plot










**Figures 7-12 (d):**

1. Use the file   plot\_energy\_loss.py
1. See lines 7 and 8 for data files location,  the data files are energy\_lost\_nocoll.dat

and energy\_lost\_coll.dat, such files are at multiple locations (obviously for different equilibrium configurations, for example

/Downloads/ISHW\_2022/data/CLASSIFICATION/QI\_Subbotin\_BETA\_0.088/POSTPR\_ENELOSS/OLD/

Check: Which code generates this file?

1. Upon running the code in Visual Studio (interactive mode), we get the following plot







**Generating the Figures of**

[Journal of Comput. Physics,](https://www.cambridge.org/core/journals/journal-of-plasma-physics) 403, 109065 (2020).

Go to folder:  /afs/itp.tugraz.at/user/khanm/Downloads/ISHW\_2019/

For orbits and cuts:


Go to your github website and look into the folder

` `[python](https://github.com/itpplasma/SIMPLE/tree/Majid/python)  >> iaea2019 >> orbits\_and\_cuts

There you will have the python code which generates various orbits (passing and banana) by choosing different values of \lambda and also starting radial coordinates. To summarize it

1. Compute orbit in cylindrical coordinates (line 153 ..)

1. Compute and plot Poincare cut (tips of banana)  in R, Z coordinates  (line 218 onward)

1. Same in topological toroidal coordinates  (line 277 onward)

**May be its possible that your interactive mode generate:**  ImportError,   

this is due to the fact that you are not in the right folder. For that you need to follow the steps

1) In interactive mode type      cd .. / .. /build    (its the place where SIMPLE can be accessed).
1) Then click somewhere on the screen of  orbits\_and\_cuts.py  and run the code in interactive mode. You will see all the desired figures.

**If you have made some changes and want to implement your changes do the following**

1) Click the source control symbol in the left most pannel
1) Click on the + sign (stage changes)
1) Type some message to let others know about your changes
1) Click on Commit
1) Lastly, click on Sync Changes


**Working with Visual Studio code on laptop using WSL**

1\. Open Visual Studio Code first

2\. Click the less/greater symbol in the bottom left corner.

3\. Connect to WSL and enter   (you can also connect to GitHub code space from here)

4\.  You should not use local uploading of any file or folder, for example to compile a code which is located in drive D of my laptop (D:\Code-SIMPLE\SIMPLE-master\SIMPLE-master\)  write the command like following

/mnt/dCode-SIMPLE/SIMPLE-master/SIMPLE-master    (Note the mnt command here)

and then follow all the steps used to compile the code

5\. At the end you can close the program easily.

6\. Remember that while plotting in my laptop some additional libraries have to be imported to have the graphical interface for the use, e.g.

I have to add the following two (in visualise3-MK)

` `14  import matplotlib   

15   matplotlib.use(‘TkAgg’)

And also had to install    sudo apt-get install python3-tk

Also the address in the basedir is written with /mnt/

Only after that it worked. So be careful while you copy the python file from Linux base system (where all such packages are already installed) to your laptop, it may not produce the graphics unless you include the required commands and installed proper packages.


SIMPLE applied to TOKAMAKS

What if I want to analyze tokamak VMEC file with SIMPLE?

It may be possible, you need to try it. First go to the website

<https://github.com/hiddenSymmetries/simsopt/tree/master/tests/test_files>

there for example you can copy the ITER file

[wout_ITERModel_reference.nc](https://github.com/hiddenSymmetries/simsopt/blob/master/tests/test_files/wout_ITERModel_reference.nc)

It has to to be an output file and the extension must be .nc.

There are also other such files in the same weblink, students can try with different comparisons there.












#
