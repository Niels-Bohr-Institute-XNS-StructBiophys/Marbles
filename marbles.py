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

has_plot = True

try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except ModuleNotFoundError:
    has_plot = False

try:
    import numpy as np
except ModuleNotFoundError:
    has_plot = False

def plot( args ):

    mpl.rcParams['figure.figsize'] = (7,8)
    mpl.rcParams['axes.labelsize']  = 16
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14
    mpl.rcParams['legend.fontsize'] = 14

    fit = np.loadtxt( f'{args.output}/best_fit.dat')
    exp = np.loadtxt( args.input )

    fig, (ax, axres) = plt.subplots(nrows = 2, sharex = False, gridspec_kw = {"height_ratios" : [5,1], "hspace" : 0} )

    ax.errorbar( exp[:,0] * 10, exp[:,1], exp[:,2], zorder = 1, fmt = 's', markersize = 6,
                  label = 'Experiment', color = '#9400d3', alpha = 0.5 )
    ax.plot( fit[:,0] * 10, fit[:,1], color = 'k', linewidth = 2, zorder = 10, label = 'Model Fit' )
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_ylabel(r'I(q) [cm$^{-1}$]')
    axres.set_xscale('log')
    axres.set_xlabel(r'q [nm$^{-1}]$')

    axres.plot( exp[:,0] * 10, (exp[:,1] - fit[:,1])/exp[:,2], color = 'k', linewidth = 2, zorder = 1 )
    axres.set_ylabel(r'$\Delta$I/$\sigma$')
    plt.tight_layout()
    plt.savefig( f'{args.output}/fit.pdf', format = 'pdf' )

def parse_inputs():

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
    parser.add_argument( '--inserted_residues', '-aa', type = int, required = True,
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
    parser.add_argument( '--schedule', '-ss', type = float, required = False, default = 0.9,
                        help = 'Simulated annealing schedule (default: 0.9)' )
    parser.add_argument( '--temperature_factor', '-t', type = float, required = False, default = 10,
                        help = 'Number by which the initial chi squared is divided to define the initial temperature (default: 10)' )

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
    parser.add_argument('--sample_info', '-si', type = str, required = False, default = None,
                        help = 'Sample information file employed in WillItFit to set the chemical properties of the nanodisc. Necessary only if the nanodisc is composed by non POPC-lipids \
                               and a belt protein different from MSP1D1')
    parser.add_argument( '--clash_distance', '-cd', type = float, required = False, default = 1.8,
                        help = 'Smaller possible distance between beads before move is rejected (default: 1.8A)' )
    # parser.add_argument( '--maximum_distance', '-md', type = float, required = False, default = 5.1,
    #                     help = 'Maximum distance between chosen beads allowed by Monte Carlo move (default: 5.1A)' )
    parser.add_argument( '--connected', '-cc', type = float, required = False, default = 5.8,
                        help = 'Maximum distance within which two beads are considered connected (default: 4.5A)' )
    parser.add_argument( '--qs_for_I0', '-qi', type = int, required = False, default = 5,
                        help = 'Number of low-q points to use to determine the value of I(0) (default: 5)' )
    parser.add_argument( '--qs_for_b', '-qb', type = int, required = False, default = 4,
                        help = 'Number of low-q points to use to determine the value of background (default: 4)' )

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    args = parse_inputs()

    # the number of loops is set as default: read the sequence file and count the residues
    if( args.loops == 0 ):
        with open( args.sequence_file ) as f:
            r = f.readlines()
        args.loops = len(r[1])

    maximum_distance = 5.1 #to be removed because it's not implemented

    if( args.convergence_acceptance == None ):
        args.convergence_acceptance = 0

    if( args.sample_info == None ):
        args.sample_info = 'None'

    cmd  = './runner ' + str(args.sequence_file) + ' ' + str(args.input) + ' ' + str(args.output) + ' ' + str(args.dmax) + \
             ' ' + str(args.passes) + ' ' + str(args.loops) + ' ' + str(args.connect_strength) + ' ' + str(args.neighbours_strength) + \
             ' ' + str(args.schedule) + ' ' + str(args.temperature_factor) + ' ' + str(args.clash_distance) + \
             ' ' + str(maximum_distance) + ' ' + str(args.connected) + ' ' + str(args.fit) + ' ' + str(args.inserted_residues) + \
             ' ' + str(args.insertion_strength) + ' ' + str(args.qs_for_I0) + ' ' + str(args.disordered_tail) + ' ' + str(args.zshift) + \
             ' ' + str(args.qs_for_b) + ' ' +  str(args.convergence_temp) + ' ' + str(args.convergence_acceptance) + \
             ' ' + str(args.intensity_stride) + ' ' + str(args.configuration_stride) + ' ' + str(args.sample_info)

    os.system( cmd )

    if has_plot:
        plot( args )
