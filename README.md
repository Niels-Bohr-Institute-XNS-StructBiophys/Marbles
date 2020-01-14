# Marbles
Author: Simone Orioli | email: simone.orioli@nbi.ku.dk    
  

C++/Python code for ab-initio shape prediction of a protein embedded in a nanodisc, given its SAXS intensity and sequence.

## Table of Contents
1. [Dependencies](#dependencies)
2. [Dependencies Installation](#dep-installation)
3. [Download and Compilation](#m-installation)
4. [Input Preparation](#input)
5. [Output](#output)
6. [Running Marbles](#running)
7. [Example Application](#example)
8. [Non-Default Samples](#non-default)
9. [Contact](#contacts)
10. [Acknowledgments](#acknow)

## Dependencies <a name="dependencies"></a>
Mandatory:
1. GSL <= 1.9
2. Python >= 2.7
3. NLOpt
4. WillItFit  
  
Additional:
1. Cmake >= 3.14 (for NLOpt and GSL installation)
2. Tkinter >= 8.6 (only for GUI)
3. PyMOL >= 2.3 (only for visualization)
4. Numpy >= 1.17 (only for plotting)
5. Matplotlib >= 3.1.1 (only for plotting)
6. Anaconda > 4.8 (only for Tkinter, Numpy and Matplotlib installation)

## Dependencies Installation <a name="dep-installation"></a>
### GSL
Marbles uses the GSL library to compute Legendre polynomials and generate pseudo-random numbers. From your console, you can download GSL 1.9 by running the following command:
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
and create and installation directory in the location you consider more appropriate (your `/home/`, for example)
```
mkdir /home/myuser/gsl-installed/
```
Mind that the name of the GSL installation directory can be whatever you prefer. Now you can configure the installation, selecting the installation directory as prefix:
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
Marbles doesn't employ WillItFit at run time, but rather WillItFit is used to fit the SAXS signal of a system of empty nanodiscs. WillItFit can be downloaded at
```
https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/Will_It_Fit.git
```
and its only dependency is the GSL library. The code doesn't require an installation, as it is compiled at run time. The use of WillItFit will not be covered here, and we refer the user to the original documentation, which can be found at
```
https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/Will_It_Fit/blob/master/WIFbashReadMe.txt
```

### PyMOL
Marbles doesn't employ PyMOL at run time, but rather it uses it to visualize the structure of the final model. PyMOL can be downloaded for free if you are an academic user or a student: 
```
https://pymol.org/2/
```
The use and installation of PyMOL will not be covered here, and we refer the user to the original documentation, which can be found at:
```
https://pymolwiki.org/index.php/Main_Page
```

### Tkinter, Numpy and Matplotlib
Marbles employs Tkinter to run a GUI, which is however not strictly necessary and can be substituted by the parser. Numpy and Matplotlib, instead, are used at the end of the program execution to plot a picture of the best fit. We suggest to install these three libraries using Anaconda (`https://www.anaconda.com/distribution/`). To install Tkinter, run in your console
```
conda install -c anaconda tk
```
Numpy and Matplotlib come already installed with Anaconda. 

## Download and Compilation <a name="m-installation"></a>
You can download the code by cloning the repository
```
https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/Marbles.git
```
Once downloaded, open `makefile` and update the variables `GSL_LIB_PATH`, `GSL_INCLUDE_PATH`, `NLOPT_LIB_PATH` and `NLOPT_INCLUDE_PATH` by specifying the paths to GSL and NLOpt. If the two libraries are not installed on your system, refer to the previous section for a quick guide on how to install them. At this point, run `make` to compile the code. The code will not be installed, and it needs to be called from within the `Marbles` directory. 

## Input Preparation <a name="input"></a>
To predict the shape of a membrane protein in a nanodisc with Marbles, you need the following information:
1. the SAXS intensity, in `cm^-1`, as a function of the scattering vector, in `A^-1`, of a solution of *empty* nanodiscs;
2. the SAXS intensity, in `cm^-1`, as a function of the scattering vector, in `A^-1`, of a solution of *loaded* nanodiscs;
3. the protein sequence;
4. the maximum diameter of the loaded nanodisc system, in `A`;
5. an estimate of the number of residues inserted in the nanodisc;
6. a `Results.wif` file, obtained from WillItFit, containing the parameters of the empty nanodisc fit.

The SAXS intensity of the empty nanodisc is used only for the determination of the nanodisc parameters through WillItFit, and will not be employed by Marbles directly. Therefore, refer to the original documentation in WillItFit for the file formatting. Instead, the SAXS intensity of the loaded nanodisc has to be provided in a file, using .dat or .txt extensions, following the format
```
q1  I1  e1
q2  I2  e2
..  ..  ..
..  ..  ..
qN  IN  eN   
```
where `qi` represent the scattering vectors, `Ii` the SAXS intensity measured at `qi` and `ei` the error on `Ii`. All comments should be removed from the file, while no particular spacing is required between columns. The protein sequence is expected to be provided in FASTA format. For example, in the case of ubiquitin the provided input file should look like
```
>|Ubiquitin(UBQ)
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
```
If some portion, or domain, of the protein is expected to be disordered and to protrude from the bottom leaflet of the nanodisc, the disordered residues should be removed from the sequence file, as Marbles will treat it on a separate footing.  

## Output<a name="output"></a>
The program outputs:
1. A `parameters.log` file with a summary of the parameters employed in the run;
2. A `penalty.dat` file saving the values of all the penalty functions at every pass during the run;
3. A `model.pdb` file, representing the best protein model generated by Marbles;
4. A `best_fit.dat` file, representing the best SAXS fit provided by Marbles;
5. An `intensity` folder, containing `*.dat` files that correspond to the computed SAXS intensities at the end of every pass;
6. A `configurations` folder, containing `*.pdb` files that correspond to the computed protein configurations at the end of every pass;
7. A `fit.pdf` file, showing the fit achieved by Marbles superimposed to the experimental SAXS intensity. 

During execution, Marbles will print on screen something similar to the following:

```

##################################
#          MARBLES v.0.1         #
#                                #
# A software for the estimation  #
# of shapes of membrane proteins #
# inserted in membrane nanodiscs #
##################################

# PRELIMINARIES
# -------------
# Results folder:          ../p450/test_nano/
# COM optimization:        X^2 = 27.0086
# Optimal sphere center:   [25, 0, 40]
# Background:              7.6e-05 (X^2 = 0.046)
# Scale factor (/1e15):    2.6

# OPTIMIZATION
# _______________________________________________________________________________________________
# Pass | Ac. R. |   Temp.   |     X     |     T     |     H     |     C     |     P     | Beads |
# -----------------------------------------------------------------------------------------------
   -1  |  ****  |     ***** |      27.0 |    6115.0 |  120175.9 |     207.3 |  126525.2 |   37  |
    0  |  0.56  |     2.701 |      77.1 |     595.0 |    2224.2 |      16.4 |    2912.6 |   16  |
    1  |  0.54  |     2.431 |      29.9 |      55.0 |     129.3 |      30.3 |     244.5 |   15  |
    2  |  0.52  |     2.188 |      19.3 |      35.0 |      61.8 |      24.9 |     141.0 |   14  |
    3  |  0.51  |     1.969 |      20.4 |      15.0 |      52.5 |      27.6 |     115.5 |   15  |

...
...
# -----------------------------------------------------------------------------------------------
```
What is reported on screen is a live diagnostic on the running Marbles process. It shows all the relevant quantities that are driving the run. The first row shows the location where outputs will be stored. The second line shows the chi squared obtained after the COM optimization, while the third line reports the position of the initial configuration at the beginning of the run. The fourth line reports the estimated backgound, while, finally, the last line shows the determined scale factor between the simulated and experimental SAXS signal. The table that follows concerns directly the optimization tun. In the first column it reports the number of passes; in the second column, it reports the acceptance ratio at the end of the pass; in the third colum, it shows the effective temperature employed throughout the ended pass; in the fourth column it shows the value of the chi squared at the end of the pass; in the fifth to eigth columns it reports the values of the penalty functions at the end of the pass, on the ninth the total penalty and on the tenth the number of beads currently inserted in the nanodisc. 

## Running Marbles <a name="running"></a>
To know all the possible options of the code, type in your console
```
python marbles.py --help
```
This command will prompt a list will all options and mandatory inputs:
```
  -h, --help            show this help message and exit
  --sequence_file SEQUENCE_FILE, -s SEQUENCE_FILE
                        Path to the protein sequence file
  --input INPUT, -i INPUT
                        Path to the file containing the experimental SAXS
                        signal. The file has to be formatted in three columns:
                        the first one for q, the second one for I(q) and the
                        last one for the error. q has to be expressed in A^-1,
                        the intensity and the error in cm^-1. The file has to
                        contain only numbers.
  --dmax DMAX, -d DMAX  Protein Dmax (in A)
  --output OUTPUT, -o OUTPUT
                        Path to directory where to store simulation results
  --fit FIT, -f FIT     Path to WillItFit output file for nanodisc best fit
  --inserted_residues INSERTED_RESIDUES, -aa INSERTED_RESIDUES
                        Number of residues to accomodate in the nanodisc
  --intensity_stride INTENSITY_STRIDE, -is INTENSITY_STRIDE
                        Stride used to print on file the calculated
                        intensities (default: 1, at the end of every pass)
  --configuration_stride CONFIGURATION_STRIDE, -cs CONFIGURATION_STRIDE
                        Stride used to print on file the calculated protein
                        configurations (default: 1, once at each pass)
  --passes PASSES, -p PASSES
                        Maximum number of Monte Carlo iterations (default:
                        100)
  --loops LOOPS, -l LOOPS
                        Number of loops per Monte Carlo iteration (default:
                        number of protein residues)
  --connect_strength CONNECT_STRENGTH, -c CONNECT_STRENGTH
                        Strength of the connectivity penalty function
                        (default: 300)
  --neighbours_strength NEIGHBOURS_STRENGTH, -n NEIGHBOURS_STRENGTH
                        Strength of the neighbours distribution penalty
                        function (default: 1)
  --insertion_strength INSERTION_STRENGTH, -ii INSERTION_STRENGTH
                        Strength of the neighbours distribution penalty
                        function (default: 5)
  --disordered_tail DISORDERED_TAIL, -dt DISORDERED_TAIL
                        Number of residues composing a disordered tail
                        protruding from the bottom leaflet of the bilayer
                        (default: 0)
  --zshift ZSHIFT, -z ZSHIFT
                        Intial shift, in A, of the protein with respect to
                        nanodisc center along the direction normal to the
                        nanodisc plane (default: 40A)
  --convergence_temp CONVERGENCE_TEMP, -ct CONVERGENCE_TEMP
                        Final temperature at which convergence is decleared
                        (default: 0.005).
  --convergence_acceptance CONVERGENCE_ACCEPTANCE, -ca CONVERGENCE_ACCEPTANCE
                        Final acceptance ratio at which convergence is
                        decleared (suggested: < 0.1). Overrides the
                        --convergence_temp and --passes option.
  --sample_info SAMPLE_INFO, -si SAMPLE_INFO
                        Sample information file employed in WillItFit to set
                        the chemical properties of the nanodisc. Necessary
                        only if the nanodisc is composed by non POPC-lipids
                        and a belt protein different from MSP1D1
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
                        connected (default: 4.5A)
  --qs_for_I0 QS_FOR_I0, -qi QS_FOR_I0
                        Number of low-q points to use to determine the value
                        of I(0) (default: 5)
  --qs_for_b QS_FOR_B, -qb QS_FOR_B
                        Number of low-q points to use to determine the value
                        of background (default: 4)
```
The user has the complete freedom to tweak all the free parameters in the code just by adjusting the inputs in the parser. However we suggest most users to rely on default parameters, as their values will work well for most of the systems. Marbles can also be ran using a GUI. To do so, type on your console:
```
python marbles_gui.py
```
This will prompt an user interface which can be filled with the same parameters shown here. 


## Example Application <a name="example"></a>
The `example/` folder contains 3 files:
1. `p450.dat`: the SAXS signal of a solution of Cytochrome P450 inserted in a nanodisc;
2. `p450.fasta.txt`: the protein sequence;
3. `Results.wif`: the best fit file obtained with WillItFit using the SAXS signal of a solution of empty nanodiscs.
It is known in the literature that the transmembrane domain of P450 is composed by 15 amino acids. As a protein Dmax, we can use approximately the half of the system's Dmax, so 60A. Given this information on the system, Marbles can be ran using the following command:
```
python marbles.py -s example/p450.fasta.txt -i example/p450.dat -d 60 -aa 15 --fit example/Results.wif -o test_run/
```
To visualize the results, it is possible to use the `plotdisc.py` script. To run it, simply type in your console
```
pymol plotdisc.py -- example/Results.wit example/test_run/model.pdb
```
The script will open PyMOL and show the protein model inserted in the nanodisc as predicted by Marbles. 

## Non-Default Samples<a name="non-default"></a>
### Presence of a Disordered Domain
As anticipated before, some samples might be composed by a disordered intracellular domain protruding from the bottom leaflet of the nanodisc. Marbles deals with these complex samples by modelling the disordered region as a random coil. To give a practical example of how the input should be provided in this case, let us assume the full sequence file of the protein of interest looks like the following:
```
>|InterestingProtein
[...]LIFAGKQLEDGRTLSDYNIQKESTLHXXXXXXXXXXXXXX
```
where `[...]` is used as a wildcard for a part of the sequence that is not shown and `XXXXXXXXXXXXXX` represents the disordered region that we want to model as a random coil. First, we remove this part of the sequence from the file:
```
>|InterestingProtein
[...]LIFAGKQLEDGRTLSDYNIQKESTLH
```
At this point, to run Marbles is necessary to specify, among the other parameters, also the length of the disordered region using the `-dt` flag. So a minimal running script would look like this:
```
python marbles.py -s prot/prot.fasta.txt -i prot/saxs.dat -d 50 -aa 20 --fit prot/Results.wif -dt 14 -o test_run/
```  

### Non-Default Nanodiscs
Marbles is programmed to interpret the input of WillItFit as coming from the fit of a POPC nanodisc with MSP1D1 belt protein. If this is not true, and the nanodisc fit has been carried out using non default parameters, it is possible to tell Marbles to override the defaults by loading the sample information file employed in WillItFit. This can be done using the `--sample_info` flag.  

### Not Saving Intermediate Configurations and Intensities
If you are planning to run multiple Marbles jobs and want to decrease the disk memory usage by not saving the intensities and structures at every pass, you can set the `-cs` and `-is` flags to increase the output stride. In particular, if you do not want any intensity and configurations printed during the execution, besides the final results, you can set the configurations and intensities stride to a number which is much bigger than the total number of passes, e.g. `-cs 1000` and `-is 1000`.
## Contacts <a name="contacts"></a>
For bug reports, issues and whatnot contact Simone Orioli at simone.orioli[at]nbi.ku.dk. If you want to engage in scientific collaborations using Marbles, contact Lise Arleth at arleth[at]nbi.ku.dk. 

## License
The code is distributed under the GPL 3.0 licence, so it is free to be modified and redistributed within the conditions of the aforementioned license.

## Aknowledgments <a name="acknow"></a>
Thanks to N. Skar-Gislinge and S. A. R. Kynde for writing part of the routines in the code, M. N. Pedersen for useful comments, C. G. Henning Hansen for testing and suggestions and Nordforsk Fund for funding.  
