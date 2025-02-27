import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import Phylo
from utils.model_window import ModalWindow
from Bio import motifs


class MotifFileDialog(ModalWindow):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(self.parent, resizable=False)
        self.title("Motif File Selector")
        self.in_file = None
        self.in_format = None

        self.supported_formats = [
            "AlignAce",
            "ClusterBuster",
            "XMS",
            "MEME",
            "MINIMAL",
            "MAST",
            "TRANSFAC",
            "pfm-four-columns",
            "pfm-four-rows",
            "pfm",
            "jaspar",
            "sites"
        ]

        # init
        self.is_closed_succesfully = False

        # Create and place labels and entry widgets
        tk.Label(self, text="Input File:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_file_entry = tk.Entry(self, width=30)
        self.input_file_entry.grid(row=0, column=1, padx=5, pady=5)

        tk.Button(self, text="Browse", command=self.browse_input_file).grid(row=0, column=2, padx=5, pady=5)

        tk.Label(self, text="Input Format:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_format_var = tk.StringVar(self)
        self.input_format_var.set(self.supported_formats[0])
        self.input_format_menu = tk.OptionMenu(self, self.input_format_var, *self.supported_formats)
        self.input_format_menu.grid(row=1, column=1, padx=5, pady=5, sticky='w')

        # Create and place the convert button
        tk.Button(self, text="Open", command=self.read_motif).grid(row=4, column=1, padx=5, pady=20)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    def browse_input_file(self):
        filename = filedialog.askopenfilename(title="Select input file")
        self.input_file_entry.delete(0, tk.END)
        self.input_file_entry.insert(0, filename)

    def read_motif(self):
        self.in_file = self.input_file_entry.get()
        self.in_format = self.input_format_var.get()

        self.is_closed_succesfully = True
        self.destroy()
