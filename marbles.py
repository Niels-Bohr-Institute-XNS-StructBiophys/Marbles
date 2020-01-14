'''*****************************************************************************
Copyright (C) 2020  Niels Bohr Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*****************************************************************************'''


import argparse
import os

parser = argparse.ArgumentParser( description = 'Marbles v. 0.1. A software for the estimation of shapes of membrane proteins inserted in membrane nanodiscs' )

# Mandatory Inputs
parser.add_argument( '--sequence_file', '-s', type = str, required = True,
                    help = 'Path to the protein sequence file' )
parser.add_argument( '--input', '-i', type = str, required = True,
                    help = 'Path to the file containing the experimental SAXS signal. The file has to be formatted in three \
                    columns: the first one for q, the second one for I(q) and the last one for the error. q has to be expressed in A^-1, \
                    the intensity and the error in cm^-1. The file has to contain only numbers.' )
parser.add_argument( '--dmax', '-d', type = float, required = True,
                    help = 'Protein Dmax (in A)' )
parser.add_argument( '--output', '-o', type = str, required = True,
                    help = 'Path to directory where to store simulation results' )
parser.add_argument( '--fit', '-f', type = str, required = True,
                    help = 'Path to WillItFit output file for nanodisc best fit' )
parser.add_argument( '--inserted_residues', '-aa', type = int, required = False,
                    help = 'Number of residues to accomodate in the nanodisc' )

#IO options
parser.add_argument( '--intensity_stride', '-is', type = int, required = False, default = 1,
                    help = 'Stride used to print on file the calculated intensities (default: 1, at the end of every pass)' )
parser.add_argument( '--configuration_stride', '-cs', type = int, required = False, default = 1,
                    help = 'Stride used to print on file the calculated protein configurations (default: 1, once at each pass)' )

# Monte Carlo options
parser.add_argument( '--passes', '-p', type = int, required = False, default = 100,
                    help = 'Maximum number of Monte Carlo iterations (default: 100)' )
parser.add_argument( '--loops', '-l', type = int, required = False, default = 0,
                    help = 'Number of loops per Monte Carlo iteration (default: number of protein residues)' )

# Penalty function options
parser.add_argument( '--connect_strength', '-c', type = float, required = False, default = 300,
                    help = 'Strength of the connectivity penalty function (default: 300)' )
parser.add_argument( '--neighbours_strength', '-n', type = float, required = False, default = 1,
                    help = 'Strength of the neighbours distribution penalty function (default: 1)' )
parser.add_argument( '--insertion_strength', '-ii', type = float, required = False, default = 5,
                    help = 'Strength of the neighbours distribution penalty function (default: 5)' )
parser.add_argument( '--disordered_tail', '-dt', type = int, required = False, default = 0,
                    help = 'Number of residues composing a disordered tail protruding from the bottom leaflet of the bilayer (default: 0)' )
parser.add_argument( '--zshift', '-z', type = float, required = False, default = 40,
                    help = 'Intial shift, in A, of the protein with respect to nanodisc center along the direction normal to the nanodisc plane (default: 40A)' )

# Convergence options
parser.add_argument( '--convergence_temp', '-ct', type = float, required = False, default = 0.005,
                    help = 'Final temperature at which convergence is decleared (default: 0.005).' )
parser.add_argument( '--convergence_acceptance', '-ca', type = float, required = False,
                    help = 'Final acceptance ratio at which convergence is decleared (suggested: < 0.1). Overrides the --convergence_temp and --passes option.' )

# Advanced options
parser.add_argument( '--schedule', '-ss', type = float, required = False, default = 0.9,
                    help = 'Simulated annealing schedule (default: 0.9)' )
parser.add_argument( '--temperature_factor', '-t', type = float, required = False, default = 10,
                    help = 'Number by which the initial chi squared is divided to define the initial temperature (default: 10)' )
parser.add_argument( '--clash_distance', '-cd', type = float, required = False, default = 1.8,
                    help = 'Smaller possible distance between beads before move is rejected (default: 1.8A)' )
parser.add_argument( '--maximum_distance', '-md', type = float, required = False, default = 5.1,
                    help = 'Maximum distance between chosen beads allowed by Monte Carlo move (default: 5.1A)' )
parser.add_argument( '--connected', '-cc', type = float, required = False, default = 5.8,
                    help = 'Maximum distance within which two beads are considered connected (default: 4.5A)' )
parser.add_argument( '--qs_for_I0', '-qi', type = int, required = False, default = 5,
                    help = 'Number of low-q points to use to determine the value of I(0) (default: 5)' )
parser.add_argument( '--qs_for_b', '-qb', type = int, required = False, default = 4,
                    help = 'Number of low-q points to use to determine the value of background (default: 4)' )

args = parser.parse_args()

# the number of loops is set as default: read the sequence file and count the residues
if( args.loops == 0 ):
    with open( args.sequence_file ) as f:
        r = f.readlines()
    args.loops = len(r[1])

if( args.convergence_acceptance == None ):
    args.convergence_acceptance = 0

cmd  = './runner ' + str(args.sequence_file) + ' ' + str(args.input) + ' ' + str(args.output) + ' ' + str(args.dmax) + \
         ' ' + str(args.passes) + ' ' + str(args.loops) + ' ' + str(args.connect_strength) + ' ' + str(args.neighbours_strength) + \
         ' ' + str(args.schedule) + ' ' + str(args.temperature_factor) + ' ' + str(args.clash_distance) + \
         ' ' + str(args.maximum_distance) + ' ' + str(args.connected) + ' ' + str(args.fit) + ' ' + str(args.inserted_residues) + \
         ' ' + str(args.insertion_strength) + ' ' + str(args.qs_for_I0) + ' ' + str(args.disordered_tail) + ' ' + str(args.zshift) + \
         ' ' + str(args.qs_for_b) + ' ' +  str(args.convergence_temp) + ' ' + str(args.convergence_acceptance) + \
         ' ' + str(args.intensity_stride) + ' ' + str(args.configuration_stride)

os.system( cmd )
