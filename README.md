# Bead Modelling
C++ code for ab-initio shape prediction of a protein given its SAXS intensity and sequence. It can be used either for solvated proteins or for proteins inserted into membrane nanodiscs. The code allows the introduction of custom penalty functions and form factors consistent with simple coarse grained models.

## Dependencies
1. GSL <= 1.9
2. Python >= 2.7
3. NLOpt

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
Thanks to Nicholas Skar Gisline and Søren Kynde for writing part of the routines in the code, Martin Nors Pedersen for useful comments and suggestions and Nordforsk Fund for funding.  
