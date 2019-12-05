# Bead Modelling
C++ code for ab-initio shape prediction of a protein given its SAXS intensity and sequence. It can be used either for solvated proteins or for proteins inserted into membrane nanodiscs. The code allows the introduction of custom penalty functions and form factors consistent with simple coarse grained models.

## Dependencies
1. GSL <= 1.9
2. Python >= 2.7
3. NLOpt

## Installation of Dependencies
### GSL
From your console, you can download GSL 1.9 by running the following command on your console
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
Mind that the name of the gsl installation directory can be whatever you prefer. Now you can configure the installation, selecting the installation directory as prefexi:
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
instead. 


## Installation
You can download the code by cloning the repository
```
git clone https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/BeadModeling.git
```
Once downloaded, open `makefile` and update variables `GSL_LIB_PATH`, `GSL_INCLUDE_PATH`, `NLOPT_LIB_PATH` and `NLOPT_INCLUDE_PATH` by specifying the path to GSL and NLOpt. If the two libraries are not installed on your system, fefer to the corrisponding documentations for more info. At this point, run `make` to compile the code.

## Run
To run the code, type in your console
```
python beads.py --help
```
This command will show all the possible options and mandatory inputs.

## Contacts
For bug reports, issues and whatnot contact simone.orioli[at]nbi.ku.dk.

## Aknowledgments
Thanks to N. Skar-Gislinge and S. A. R. Kynde for writing part of the routines in the code, M. N. Pedersen for useful comments and suggestions and Nordforsk Fund for funding.  
