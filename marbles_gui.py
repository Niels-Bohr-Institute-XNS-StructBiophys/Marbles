try:
    import tkinter
    from tkinter import ttk, filedialog, messagebox
except ModuleNotFoundError:
    print("\n# tkinter is necessary to run Marbles GUI. You can install it via: \n# `conda install -c anaconda tk`\n")
    exit(-1)

import re
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

class GUI( ttk.Frame ):

    def __init__( self, parent, *args, **kwargs ):
        ttk.Frame.__init__( self, parent, *args, **kwargs )
        self.root = parent
        self.init_gui()

    def on_quit(self):
        """Exits program."""
        quit()

    def saxs_upload( self, event = None ):
        filename = filedialog.askopenfilename()
        self.saxs.insert( tkinter.END, filename )

    def sequence_upload( self, event = None ):
        filename = filedialog.askopenfilename()
        self.sequence.insert( tkinter.END, filename )

        seq = self.sequence.get()

        with open( seq, 'r' ) as f:
            r = f.readlines()
        self.nloops.insert( tkinter.END, str(len(r[1])) )

    def fit_upload( self, event = None ):
        filename = filedialog.askopenfilename()
        self.fit.insert( tkinter.END, filename )

    def sinfo_upload( self, event = None ):
        filename = filedialog.askopenfilename()
        self.sampinf.insert( tkinter.END, filename )

    def set_out( self, event = None ):
        filename = filedialog.askdirectory()
        self.out.insert( tkinter.END, filename )

    def error( self, string, message ):
        if( string == '' and not self.error_called ):
            self.can_run      = False
            self.error_called = True
            messagebox.showerror("Error", message )
        else:
            pass

    def plot( self ):

        mpl.rcParams['figure.figsize'] = (7,8)
        mpl.rcParams['axes.labelsize']  = 16
        mpl.rcParams['xtick.labelsize'] = 14
        mpl.rcParams['ytick.labelsize'] = 14
        mpl.rcParams['legend.fontsize'] = 14

        fit = np.loadtxt( f'{self.output}/best_fit.dat')
        exp = np.loadtxt( self.input )

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
        plt.savefig( f'{self.output}/fit.pdf', format = 'pdf' )

    def run(self):

        self.can_run           = True
        self.error_called      = False
        self.input             = self.saxs.get()
        self.output            = self.out.get()
        sequence_file          = self.sequence.get()
        dmax                   = self.dmax.get()
        fit                    = self.fit.get()
        inserted_residues      = self.insert.get()

        self.error( input,             "SAXS intensity file is mandatory." )
        self.error( sequence_file,     "Sequence file is mandatory." )
        self.error( fit,               "WillItFit file is mandatory." )
        self.error( output,            "Output directory is mandatory." )
        self.error( dmax,              "Dmax is mandatory." )
        self.error( inserted_residues, "Number of inserted residues is mandatory." )

        intensity_stride       = self.istride.get()
        configuration_stride   = self.cstride.get()

        self.error( intensity_stride,     "Intensity stride field is empty. Default value is 1." )
        self.error( configuration_stride, "Configuration stride field is empty. Default value is 1." )

        passes                 = self.npasses.get()
        loops                  = self.nloops.get()
        schedule               = self.cooling.get()
        temperature_factor     = self.tfact.get()

        self.error( passes,             "Number of passes field is empty. Default value is 99." )
        self.error( loops,              "Loops per pass field is empty. Default value is the number of protein residues." )
        self.error( schedule,           "Cooling schedule field is empty. Default value is 0.9.")
        self.error( temperature_factor, "Temperature factor field is empty. Default value is 10.")

        connect_strength       = self.c0.get()
        neighbours_strength    = self.h0.get()
        insertion_strength     = self.t0.get()
        disordered_tail        = self.disord.get()
        zshift                 = self.shift.get()

        self.error( connect_strength,    "c0 field is empty. Default value is 99." )
        self.error( neighbours_strength, "h0 field is empty. Default value is the number of protein residues." )
        self.error( insertion_strength,  "t0 field is empty. Default value is 0.9.")
        self.error( disordered_tail,     "Disordered residues field is empty. Default value is 0.")
        self.error( zshift,              "Z shift field is empty. Default value is 40 A.")

        convergence_temp       = self.beta.get()
        convergence_acceptance = self.accr.get()

        self.error( convergence_temp, "Beta field is empty. Default value is 0.005." )

        sample_info            = self.sampinf.get()
        clash_distance         = self.cldist.get()
        connected              = self.codist.get()
        qs_for_I0              = self.qI0.get()
        qs_for_b               = self.qB.get()

        self.error( connected,          "Connection field is empty. Default value is 5.8 A." )
        self.error( clash_distance,     "Clash distance field is empty. Default value is 1.8 A." )
        self.error( qs_for_I0,          "q points for I(0) field is empty. Default value is 5.")
        self.error( qs_for_b,           "q points for B field is empty. Default value is 4.")

        if( convergence_acceptance == "" ):
            convergence_acceptance = 0

        if( sample_info == "" ):
            sample_info = 'None'

        maximum_distance = 5.1

        if self.can_run:
            cmd  = './runner ' + str(sequence_file)      + ' ' + str(input)                  + ' ' + str(output)            + ' ' + str(dmax) +                 \
                           ' ' + str(passes)             + ' ' + str(loops)                  + ' ' + str(connect_strength)  + ' ' + str(neighbours_strength) +  \
                           ' ' + str(schedule)           + ' ' + str(temperature_factor)     + ' ' + str(clash_distance)    + ' ' + str(maximum_distance)   +   \
                           ' ' + str(connected)          + ' ' + str(fit)                    + ' ' + str(inserted_residues) + ' ' + str(insertion_strength) +   \
                           ' ' + str(qs_for_I0)          + ' ' + str(disordered_tail)        + ' ' + str(zshift)            + ' ' + str(qs_for_b) +             \
                           ' ' + str(convergence_temp)   + ' ' + str(convergence_acceptance) + ' ' + str(intensity_stride)  + ' ' + str(configuration_stride) + \
                           ' ' + str(sample_info)

            os.system( cmd )

            if has_plot:
                plot()

    def init_gui(self):

        """Builds GUI."""

        self.root.title( 'Marbles v0.1: Ab Initio Protein Shape Prediction' )
        self.root.option_add( '*tearOff', 'FALSE' )

        self.grid( column = 0, row = 0, sticky = 'nsew' )

        self.menubar = tkinter.Menu( self.root )

        self.menu_file = tkinter.Menu( self.menubar )
        self.menu_file.add_command( label = 'Exit', command = self.on_quit )

        self.menu_edit = tkinter.Menu( self.menubar )

        self.menubar.add_cascade( menu = self.menu_file, label = 'File' )
        self.menubar.add_cascade( menu = self.menu_edit, label = 'Edit' )

        self.root.config( menu = self.menubar )

        self.saxs  = ttk.Entry(  self, width = 30 )
        self.bsaxs = ttk.Button( self, text = 'Browse', command = self.saxs_upload )
        self.saxs.grid(  column = 1, row = 2, columnspan = 2 )
        self.bsaxs.grid( column = 3, row = 2, columnspan = 1 )

        self.sequence  = ttk.Entry(self, width = 30 )
        self.bsequence = ttk.Button(self, text = 'Browse', command = self.sequence_upload )
        self.sequence.grid(  column = 1, row = 3, columnspan = 2 )
        self.bsequence.grid( column = 3, row = 3, columnspan = 1 )

        self.fit  = ttk.Entry(  self, width = 30 )
        self.bfit = ttk.Button( self, text = 'Browse', command = self.fit_upload )
        self.fit.grid(  column = 1, row = 4, columnspan = 2 )
        self.bfit.grid( column = 3, row = 4, columnspan = 1 )

        self.out  = ttk.Entry(  self, width = 30 )
        self.bout = ttk.Button( self, text = 'Browse', command = self.set_out )
        self.out.grid(  column = 1, row = 5, columnspan = 2 )
        self.bout.grid( column = 3, row = 5, columnspan = 1 )

        self.dmax = ttk.Entry( self, width = 10 )
        self.dmax.grid( column = 1, row = 6 )
        self.insert = ttk.Entry( self, width = 10 )
        self.insert.grid( column = 3, row = 6 )
        #-----------------------------------------------------------------------

        self.istride = ttk.Entry( self, width = 10 )
        self.istride.grid( column = 1, row = 10 )
        self.istride.insert( tkinter.END, '1' )
        self.cstride = ttk.Entry( self, width = 10 )
        self.cstride.grid( column = 3, row = 10 )
        self.cstride.insert( tkinter.END, '1' )
        #-----------------------------------------------------------------------

        self.npasses = ttk.Entry( self, width = 10, text = '3' )
        self.npasses.grid( column = 1, row = 14 )
        self.npasses.insert( tkinter.END, '99' )
        self.nloops = ttk.Entry( self, width = 10 )
        self.nloops.grid( column = 3, row = 14 )

        self.cooling = ttk.Entry( self, width = 10 )
        self.cooling.grid( column = 1, row = 15 )
        self.cooling.insert( tkinter.END, '0.9' )
        self.tfact = ttk.Entry( self, width = 10 )
        self.tfact.grid( column = 3, row = 15 )
        self.tfact.insert( tkinter.END, '10' )
        #-----------------------------------------------------------------------

        self.c0 = ttk.Entry( self, width = 10 )
        self.c0.grid( column = 1, row = 19 )
        self.c0.insert( tkinter.END, '300' )
        self.h0 = ttk.Entry( self, width = 10 )
        self.h0.grid( column = 3, row = 19 )
        self.h0.insert( tkinter.END, '1' )

        self.t0 = ttk.Entry( self, width = 10 )
        self.t0.grid( column = 1, row = 20 )
        self.t0.insert( tkinter.END, '5' )

        self.disord = ttk.Entry( self, width = 10 )
        self.disord.grid( column = 3, row = 20 )
        self.disord.insert( tkinter.END, '0' )

        self.shift = ttk.Entry( self, width = 10 )
        self.shift.grid( column = 1, row = 21 )
        self.shift.insert( tkinter.END, '40' )
        #-----------------------------------------------------------------------

        self.beta = ttk.Entry( self, width = 10 )
        self.beta.grid( column = 1, row = 25 )
        self.beta.insert( tkinter.END, '0.005' )
        self.accr = ttk.Entry( self, width = 10 )
        self.accr.grid( column = 3, row = 25 )
        #-----------------------------------------------------------------------

        self.sampinf  = ttk.Entry(self, width = 30 )
        self.bsampinf = ttk.Button(self, text = 'Browse', command = self.sinfo_upload )
        self.sampinf.grid(  column = 1, row = 29, columnspan = 2 )
        self.bsampinf.grid( column = 3, row = 29, columnspan = 1 )

        self.cldist = ttk.Entry( self, width = 10 )
        self.cldist.grid( column = 1, row = 30 )
        self.cldist.insert( tkinter.END, '1.8' )
        self.codist = ttk.Entry( self, width = 10 )
        self.codist.grid( column = 3, row = 30 )
        self.codist.insert( tkinter.END, '5.8' )

        self.qI0 = ttk.Entry( self, width = 10 )
        self.qI0.grid( column = 1, row = 31 )
        self.qI0.insert( tkinter.END, '5' )
        self.qB = ttk.Entry( self, width = 10 )
        self.qB.grid( column = 3, row = 31 )
        self.qB.insert( tkinter.END, '4' )
        #-----------------------------------------------------------------------

        self.calc_button = ttk.Button( self, text = 'Run', command = self.run )
        self.calc_button.grid( column = 0, row = 33, columnspan = 4 )

        # Labels that remain constant throughout execution.
        ttk.Label( self, text = 'Mandatory Options'       ).grid( column = 0, row = 0, columnspan = 4 )
        ttk.Separator( self, orient = 'horizontal'        ).grid( column = 0, row = 1, columnspan = 4, sticky = 'ew' )

        ttk.Label( self, text = 'SAXS Intensity'          ).grid( column = 0, row = 2, sticky = 'w' )
        ttk.Label( self, text = 'Sequence file'           ).grid( column = 0, row = 3, sticky = 'w' )
        ttk.Label( self, text = 'WillItFit file'          ).grid( column = 0, row = 4, sticky = 'w' )
        ttk.Label( self, text = 'Output'                  ).grid( column = 0, row = 5, sticky = 'w' )
        ttk.Label( self, text = 'Dmax'                    ).grid( column = 0, row = 6, sticky = 'w' )
        ttk.Label( self, text = 'Inserted residues'       ).grid( column = 2, row = 6, sticky = 'w')
        #-----------------------------------------------------------------------

        ttk.Label( self, text = ''                        ).grid( column = 0, row = 7, columnspan = 4 )
        ttk.Label( self, text = 'IO Options'              ).grid( column = 0, row = 8, columnspan = 4 )
        ttk.Separator( self, orient = 'horizontal'        ).grid( column = 0, row = 9, columnspan = 4, sticky = 'ew' )

        ttk.Label( self, text = 'Intensity stride'        ).grid( column = 0, row = 10, sticky = 'w' )
        ttk.Label( self, text = 'Configuration stride'    ).grid( column = 2, row = 10, sticky = 'w' )
        #-----------------------------------------------------------------------

        ttk.Label( self, text = ''                        ).grid( column = 0, row = 11, columnspan = 4 )
        ttk.Label( self, text = 'Monte Carlo Options'     ).grid( column = 0, row = 12, columnspan = 4 )
        ttk.Separator( self, orient = 'horizontal'        ).grid( column = 0, row = 13, columnspan = 4, sticky = 'ew' )

        ttk.Label( self, text = 'Number of passes'        ).grid( column = 0, row = 14, sticky = 'w' )
        ttk.Label( self, text = 'Loops per pass'          ).grid( column = 2, row = 14, sticky = 'w' )
        ttk.Label( self, text = 'Cooling schedule'        ).grid( column = 0, row = 15, sticky = 'w' )
        ttk.Label( self, text = 'Temperature factor'      ).grid( column = 2, row = 15, sticky = 'w' )
        #-----------------------------------------------------------------------

        ttk.Label( self, text = ''                        ).grid( column = 0, row = 16, columnspan = 4 )
        ttk.Label( self, text = 'Penalty Options'         ).grid( column = 0, row = 17, columnspan = 4 )
        ttk.Separator( self, orient = 'horizontal'        ).grid( column = 0, row = 18, columnspan = 4, sticky = 'ew' )

        ttk.Label( self, text = 'c0'                      ).grid( column = 0, row = 19, sticky = 'w' )
        ttk.Label( self, text = 'h0'                      ).grid( column = 2, row = 19, sticky = 'w' )
        ttk.Label( self, text = 't0'                      ).grid( column = 0, row = 20, sticky = 'w' )
        ttk.Label( self, text = 'Disordered residues'     ).grid( column = 2, row = 20, sticky = 'w' )
        ttk.Label( self, text = 'Z Shift (A)'             ).grid( column = 0, row = 21, sticky = 'w' )
        #-----------------------------------------------------------------------

        ttk.Label( self, text = ''                        ).grid( column = 0, row = 22, columnspan = 4 )
        ttk.Label( self, text = 'Convergence Options'     ).grid( column = 0, row = 23, columnspan = 4 )
        ttk.Separator( self, orient = 'horizontal'        ).grid( column = 0, row = 24, columnspan = 4, sticky = 'ew' )

        ttk.Label( self, text = 'Beta'                    ).grid( column = 0, row = 25, sticky = 'w' )
        ttk.Label( self, text = 'Acceptance ratio'        ).grid( column = 2, row = 25, sticky = 'w' )
        #-----------------------------------------------------------------------

        ttk.Label( self, text = ''                        ).grid( column = 0, row = 26, columnspan = 4 )
        ttk.Label( self, text = 'Advanced Options'        ).grid( column = 0, row = 27, columnspan = 4 )
        ttk.Separator( self, orient = 'horizontal'        ).grid( column = 0, row = 28, columnspan = 4, sticky = 'ew' )

        ttk.Label( self, text = 'Sample information'      ).grid( column = 0, row = 29, sticky = 'w' )
        ttk.Label( self, text = 'Clash distance (A)'      ).grid( column = 2, row = 30, sticky = 'w' )
        ttk.Label( self, text = 'Connection distance (A)' ).grid( column = 0, row = 30, sticky = 'w' )
        ttk.Label( self, text = 'q points for I(0)'       ).grid( column = 0, row = 31, sticky = 'w' )
        ttk.Label( self, text = 'q points for B'          ).grid( column = 2, row = 31, sticky = 'w' )
        #-----------------------------------------------------------------------


        for child in self.winfo_children():
            child.grid_configure( padx = 5, pady = 5 )

if __name__ == '__main__':
    root = tkinter.Tk()
    GUI( root )
    root.mainloop()
