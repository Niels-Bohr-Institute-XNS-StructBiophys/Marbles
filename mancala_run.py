import tkinter
from tkinter import ttk
from tkinter import filedialog
import re

class Adder(ttk.Frame):
    """The adders gui and functions."""
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.root = parent
        self.init_gui()

    def on_quit(self):
        """Exits program."""
        quit()

    def UploadAction(self, event=None):
        filename = filedialog.askopenfilename()
        print('Selected:', filename)

    def calculate(self):
        """Calculates the sum of the two inputted numbers."""
        #num1 = int(self.num1_entry.get())
        #num2 = int(self.num2_entry.get())
        #num3 = num1 + num2
        #self.answer_label['text'] = num3

        num1 = int(self.num1_entry.get())
        num2 = int(self.num2_entry.get())
        num3 = int(self.num3_entry.get())
        num4 = int(self.num4_entry.get())

        print("HEI YO! Complex operations here: ", num1, num2, num3, num4 )

    def init_gui(self):
        """Builds GUI."""
        self.root.title('Mancala: SAXS Data Bead Modeling Routine')
        self.root.option_add('*tearOff', 'FALSE')

        self.grid(column=0, row=0, sticky='nsew')

        self.menubar = tkinter.Menu(self.root)

        self.menu_file = tkinter.Menu(self.menubar)
        self.menu_file.add_command(label='Exit', command=self.on_quit)

        self.menu_edit = tkinter.Menu(self.menubar)

        self.menubar.add_cascade(menu=self.menu_file, label='File')
        self.menubar.add_cascade(menu=self.menu_edit, label='Edit')

        self.root.config(menu=self.menubar)

        self.data_button = ttk.Button(self, text='Browse',command=self.UploadAction)
        self.data_button.grid(column=1, row=2, columnspan=1)

        # button = tkinter.Button(root, text='Open', command=self.UploadAction, width=5)
        # button.grid(column=1, row=2)

        #self.datafile = ttk.Entry(self, width=5)
        #self.datafile.grid(column=1, row = 2)

        #self.fit = ttk.Entry(self, width=5)
        #self.fit.grid(column=3, row=2)

        self.fit_button = ttk.Button(self, text='Browse',command=self.UploadAction)
        self.fit_button.grid(column=3, row=2, columnspan=1)

        #self.sequence = ttk.Entry(self, width=5)
        #self.sequence.grid(column=1, row = 3)

        self.sequence_button = ttk.Button(self, text='Browse',command=self.UploadAction)
        self.sequence_button.grid(column=1, row=3, columnspan=1)

        #self.pdb = ttk.Entry(self, width=5)
        #self.pdb.grid(column=3, row=3)

        self.pdb_button = ttk.Button(self, text='Browse',command=self.UploadAction)
        self.pdb_button.grid(column=3, row=3, columnspan=1)

        self.dmax = ttk.Entry(self, width=8)
        self.dmax.grid(column=1, row=4)

        self.npasses = ttk.Entry(self, width=8)
        self.npasses.insert(tkinter.END, '80')
        self.npasses.grid(column=3, row=4)

        self.loops_per_pass = ttk.Entry(self, width=8)
        self.loops_per_pass.insert(tkinter.END, 'def')
        self.loops_per_pass.grid(column=1, row=5)

        self.ll = ttk.Entry(self, width=8)
        self.ll.insert(tkinter.END, '100')
        self.ll.grid(column=3, row=5)

        self.connect = ttk.Entry(self, width=8)
        self.connect.insert(tkinter.END, '30000')
        self.connect.grid(column=1, row=6)

        self.verbose = ttk.Entry(self, width=8)
        self.verbose.insert(tkinter.END, 'yes')
        self.verbose.grid(column=3, row=6)



        self.calc_button = ttk.Button(self, text='Run',
                command=self.calculate)
        self.calc_button.grid(column=0, row=10, columnspan=4)

        #self.answer_frame = ttk.LabelFrame(self, text='Answer',
        #        height=100)
        #self.answer_frame.grid(column=0, row=4, columnspan=4, sticky='nesw')

        #self.answer_label = ttk.Label(self.answer_frame, text='')
        #self.answer_label.grid(column=0, row=0)

        # Labels that remain constant throughout execution.
        ttk.Label(self, text='Input parameters').grid(column=0, row=0,columnspan=4)
        ttk.Separator(self, orient='horizontal').grid(column=0,row=1, columnspan=4, sticky='ew')
        ttk.Label(self, text='Data file (.rad)').grid(column=0, row=2,sticky='w')
        ttk.Label(self, text='Fit file (*_bf.dat)').grid(column=2, row=2,sticky='w')
        ttk.Label(self, text='Sequence file').grid(column=0, row=3,sticky='w')
        ttk.Label(self, text='Sequence from PDB').grid(column=2, row=3,sticky='w')
        ttk.Label(self, text='Dmax').grid(column=0, row=4,sticky='w')
        ttk.Label(self, text='Number of passes').grid(column=2, row=4,sticky='w')
        ttk.Label(self, text='Loops per pass').grid(column=0, row=5,sticky='w')
        ttk.Label(self, text='Lambda').grid(column=2, row=5,sticky='w')
        ttk.Label(self, text='Connect').grid(column=0, row=6,sticky='w')
        ttk.Label(self, text='Verbose').grid(column=2, row=6,sticky='w')

        ttk.Separator(self, orient='horizontal').grid(column=0,row=7, columnspan=4, sticky='ew')
        ttk.Label(self, text='Advanced parameters').grid(column=0, row=8,columnspan=4)

        for child in self.winfo_children():
            child.grid_configure(padx=5, pady=5)

if __name__ == '__main__':
    root = tkinter.Tk()
    Adder(root)
    root.mainloop()
