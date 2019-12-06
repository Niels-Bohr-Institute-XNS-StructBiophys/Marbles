# Bead Modelling
C++ code for ab-initio shape prediction of a protein given its SAXS intensity and sequence. It can be used either for solvated proteins or for proteins inserted into membrane nanodiscs. The code allows the introduction of custom penalty functions and form factors consistent with simple coarse grained models.

## Dependencies
1. GSL <= 1.9
2. Python >= 2.7
3. Cmake
4. NLOpt
5. WillItFit

## Installation of Dependencies
### GSL
Marbles uses the GSL library to compute legendre polynomials and random numbers. From your console, you can download GSL 1.9 by running the following command:
```
wget https://ftp.gnu.org/gnu/gsl/gsl-1.9.tar.gz
```
Place the GSL folder in your home directory and unpack it by running
```
tar xzf gsl-1.9.tar.gz
```
Change to this directory
```
cd gsl-1.9
```
and create and installation directory in the location you consider more appropriate (your /home/, for example)
```
mkdir /home/myuser/gsl-installed/
```
Mind that the name of the gsl installation directory can be whatever you prefer. Now you can configure the installation, selecting the installation directory as prefix:
```
./configure --prefix=/home/myuser/gsl-installed/
```
At this point you can compile the library by running
```
make
```
and check that the installation went fine with
```
make check
```
You are now ready for installation
```
make install
```
Beware that if you didn't specify any prefix during configuration, you will need super user privileges to install the library. If so, just run
```
sudo make install
```
instead. This is the list of the minimal steps required to make Marbles work, but we refer the user to the official GSL page, `https://www.gnu.org/software/gsl/` for further information about the library. 

### NLOpt
Marbles employes NLOpt library to run some simple fits. From your console, you can download NLOpt running the following command:
```
wget https://github.com/stevengj/nlopt/archive/v2.6.1.tar.gz
```
Place the NLOpt folder in your home directory and unpack it by running
```
tar xzf v2.6.1.tar.gz
```
Change to this directory
```
cd nlopt-2.6.1
```
Create a build directory
```
mkdir build
```
and change to it
```
cd build
```
At this point, configure the installation using cmake, specifying `/home/myuser/nlopt-installed/` as the installation directory
```
cmake .. -DCMAKE_INSTALL_PREFIX=/home/myuser/nlopt-installed/
```
Now compile the library
```
make
```
and install it
```
make install
```
If no prefix has been specified at configuration level, a default directory will be chosen and you will need super user privileges to install it. In that case, run 
```
sudo make install
```
instead. This is the list of the minimal steps required to make Marbles work, but we refer the user to the official NLOpt page, `https://nlopt.readthedocs.io/en/latest/` for further information about the library. 
### WillItFit
Marbles doesn't employ WillItFit at run time, but rather WillItFit is used to fit the SAXS signal of the bare nanodisc. WillItFit can be downloaded at
```
https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/Will_It_Fit.git
```
and its only dependency is the GSL library. The use of WillItFit will not be covered here, and we refer the user to the original documentation, which can be found at
```
https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/Will_It_Fit/blob/master/WIFbashReadMe.txt
```

## Marbles Download and Compilation
You can download the code by cloning the repository
```
git clone https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/BeadModeling.git
```
Once downloaded, open `makefile` and update the variables `GSL_LIB_PATH`, `GSL_INCLUDE_PATH`, `NLOPT_LIB_PATH` and `NLOPT_INCLUDE_PATH` by specifying the paths to GSL and NLOpt. If the two libraries are not installed on your system, refer to the previous section for a quick guide on how to install them. At this point, run `make` to compile the code. The code will not be installed, and it needs to be called from the `marbles` directory. 

## Input Preparation
To predict the shape of a membrane protein in a nanodisc with Marbles, you need the following information:
1. the SAXS intensity, in cm^-1, as a function of the scattering vector, in A^-1, of a solution of empty nanodiscs;
2. the SAXS intensity, in cm^-1, as a function of the scattering vector, in A^-1, of a solution of loaded nanodiscs;
3. the protein sequence;
4. the Dmax of the loaded nanodisc system;
5. an estimate of the number of residues inserted in the nanodisc;
6. a Results.wif file, obtained from WillItFit, containing the parameters of the empty nanodisc fit.
The SAXS intensity of the empty nanodisc is used only for the determination of the nanodisc parameters through WillItFit, and will not be employed by Marbles directly. Therefore, refer to the original documentation in WillItFit for the file formatting. Instead, the SAXS intensity of the loaded nanodisc has to be provided in a file, using .dat or .txt extensions, following the format
```
q1  I1  e1
q2  I2  e2
..  ..  ..
..  ..  ..
qN  IN  eN   
```
where qi represents the scattering vector, Ii the SAXS intensity measured at qi and ei the error on Ii. All comments should be removed from the file, while no particular spacing is required between columns. The protein sequence is expected to be provided in FASTA format. If some portion, or domain, of the protein is expected to be disordered and to protrude from the bottom leaflet of the nanodisc, we suggest the user to remove this part of the sequence from the FASTA file, as Marbles will treat it on a separate footing.  

## Running Marbles
To know all the possible options of the code, type in your console
```
python marbles.py --help
```
This command will prompt a list will all options and mandatory inputs:
```
  -h, --help            show this help message and exit
  --sequence_file SEQUENCE_FILE, -s SEQUENCE_FILE
                        Path to the protein sequence file
  --with_nanodisc, -w   Use to run the simulation in the presence of a
                        nanodisc; do not use for protein-only run.
  --input INPUT, -i INPUT
                        Path to the file containing the experimental SAXS
                        signal. The file has to be formatted in three columns:
                        the first one for q, the second one for I(q) and the
                        last one for the error. q has to be expressed in A^-1,
                        the intensity and the error in cm^-1. The file has to
                        contain only numbers.
  --dmax DMAX, -d DMAX  Protein Dmax (in A)
  --loops LOOPS, -l LOOPS
                        Number of loops per Monte Carlo iteration
  --output OUTPUT, -o OUTPUT
                        Path to directory where to store simulation results
  --fit FIT, -f FIT     Path to WillItFit output file for nanodisc best fit
  --passes PASSES, -p PASSES
                        Number of Monte Carlo iterations (default: 100)
  --connect_strength CONNECT_STRENGTH, -c CONNECT_STRENGTH
                        Strength of the connectivity penalty function
                        (default: 300)
  --neighbours_strength NEIGHBOURS_STRENGTH, -n NEIGHBOURS_STRENGTH
                        Strength of the neighbours distribution penalty
                        function (default: 1)
  --inserted_residues INSERTED_RESIDUES, -aa INSERTED_RESIDUES
                        Number of residues to accomodate in the nanodisc
  --insertion_strength INSERTION_STRENGTH, -ii INSERTION_STRENGTH
                        Strength of the neighbours distribution penalty
                        function (default: 5)
  --schedule SCHEDULE, -ss SCHEDULE
                        Simulated annealing schedule (default: 0.9)
  --temperature_factor TEMPERATURE_FACTOR, -t TEMPERATURE_FACTOR
                        Number by which the initial chi squared is divided to
                        define the initial temperature (default: 10)
  --clash_distance CLASH_DISTANCE, -cd CLASH_DISTANCE
                        Smaller possible distance between beads before move is
                        rejected (default: 1.8A)
  --maximum_distance MAXIMUM_DISTANCE, -md MAXIMUM_DISTANCE
                        Maximum distance between chosen beads allowed by Monte
                        Carlo move (default: 5.1A)
  --connected CONNECTED, -cc CONNECTED
                        Maximum distance within which two beads are considered
                        connected (default: 5.81A)
  --qs_for_I0 QS_FOR_I0, -qi QS_FOR_I0
                        Number of low-q points to use to determine the value
                        of I(0) (default: 5)
```
The user has the complete freedom to tweak all the free parameters in the code just by adjusting the inputs in the parser. However, inexperienced users can ignore most of the flags, as the default values will work well for most of the systems. A minimal running command for a protein in solution composed by 300 amino acids and with a Dmax of 50A would look like this:
```
python marbles.py -s /my/input/prot.fasta.txt -i /my/input/saxs.dat -d 50 -l 300 -o /my/output/marbles_run/
```
The `/my/output/` should exist already and will not be created. Instead, a minimal running command for a protein in a nanodisc composed by 300 amino acids and with a Dmax of 50A, with 15 amino acids inserted in the bilayer, would look like this:
```
python marbles.py -w -s /my/input/prot.fasta.txt -i /my/input/saxs.dat -d 50 -l 300 -aa 15 --fit /my/input/Results.wif -o /my/output/marbles_run/
```

## Contacts
For bug reports, issues and whatnot contact simone.orioli[at]nbi.ku.dk.

## Aknowledgments
Thanks to N. Skar-Gislinge and S. A. R. Kynde for writing part of the routines in the code, M. N. Pedersen for useful comments and suggestions and Nordforsk Fund for funding.  
