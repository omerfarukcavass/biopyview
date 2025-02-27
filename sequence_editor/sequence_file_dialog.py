import os
import tkinter as tk
from tkinter import filedialog, messagebox
from utils.model_window import ModalWindow
from Bio import SeqIO
from tkinter import ttk

from utils.utils import is_valid_path
from utils.windows import TKMessageBox


class SequenceOpenFileDialog(ModalWindow):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(self.parent, resizable=False)
        self.title("Sequence File Dialog")

        self.supported_formats = list(SeqIO._FormatToIterator.keys())  # All return a SeqRecord iterator.

        # self.supported_formats = [
        #     "abi",
        #     "abi-trim",
        #     "ace",
        #     "fasta",
        #     "fasta-2line",
        #     "ig",
        #     "embl",
        #     "embl-cds",
        #     "gb",
        #     "gck",  # alias for genbank
        #     "genbank",
        #     "genbank-cds",
        #     "imgt",
        #     "nib",
        #     "cif-seqres",
        #     "cif-atom",
        #     "pdb-atom",
        #     "pdb-seqres",
        #     "phd",
        #     "pir",
        #     "fastq",
        #     "fastq-sanger",  # alias for fastq
        #     "fastq-solexa",
        #     "fastq-illumina",
        #     "qual",
        #     "seqxml",
        #     "sff",
        #     "snapgene",
        #     "sff-trim",
        #     "swiss",
        #     "tab",
        #     "twobit",
        #     "uniprot-xml",
        #     "xdna"
        # ]

        self.zip_formats = ["-","gz", "bz2"]

        # init
        self.is_closed_succesfully = False
        self.in_file = None
        self.in_format = None
        self.zip_format = None

        # Create and place labels and entry widgets
        tk.Label(self, text="Input File:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_file_entry = tk.Entry(self, width=30)
        self.input_file_entry.grid(row=0, column=1, padx=5, pady=5, sticky='we')

        tk.Button(self, text="Browse", command=self.browse_input_file).grid(row=0, column=4, padx=5, pady=5)

        tk.Label(self, text="Input Format:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_format_var = tk.StringVar(self)
        self.input_format_var.set("fasta")
        self.input_format_menu = tk.OptionMenu(self, self.input_format_var, *self.supported_formats)
        self.input_format_menu.grid(row=1, column=1, padx=5, pady=5,  sticky='w')

        tk.Label(self, text="Zip Format:").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.zip_format_var = tk.StringVar(self)
        self.zip_format_var.set(self.zip_formats[0])
        self.zip_format_menu = tk.OptionMenu(self, self.zip_format_var, *self.zip_formats)
        self.zip_format_menu.grid(row=2, column=1, padx=5, pady=5, sticky='w')

        # Create and place the convert button
        tk.Button(self, text="Open", command=self.open_file).grid(row=3, column=1, padx=5, pady=20)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    def disable_menu(self):
        pass

    def enable_menu(self):
        pass

    def browse_input_file(self):
        filename = filedialog.askopenfilename(title="Select input file")
        self.input_file_entry.delete(0, tk.END)
        self.input_file_entry.insert(0, filename)

    def open_file(self):
        self.in_file = self.input_file_entry.get()
        self.in_format = self.input_format_var.get()
        self.zip_format = self.zip_format_var.get()

        if not is_valid_path(self.in_file):
            TKMessageBox("Error", f"Please select a valid path for input file.", messagebox.showerror, self).show()
            return

        self.is_closed_succesfully = True
        self.destroy()


class SequenceOpenIndexFileDialog(ModalWindow):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(self.parent, resizable=False)
        self.title("Sequence Index File Dialog")

        from Bio.SeqIO._index import _FormatToRandomAccess
        self.supported_formats = list(_FormatToRandomAccess.keys())

        # self.supported_formats = [
        #     "ace",
        #     "embl",
        #     "fasta",
        #     "fastq",
        #     "fastq-sanger",
        #     "fastq-solexa",
        #     "fastq-illumina",
        #     "genbank",
        #     "gb",
        #     "ig",
        #     "imgt",
        #     "phd",
        #     "pir",
        #     "sff",
        #     "sff-trim",
        #     "swiss",
        #     "tab",
        #     "qual",
        #     "uniprot-xml"
        # ]

        # init
        self.is_closed_succesfully = False
        self.in_file = None
        self.in_format = None

        # Create and place labels and entry widgets
        tk.Label(self, text="Input File:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_file_entry = tk.Entry(self, width=30)
        self.input_file_entry.grid(row=0, column=1, padx=5, pady=5, sticky='we')

        tk.Button(self, text="Browse", command=self.browse_input_file).grid(row=0, column=4, padx=5, pady=5)

        tk.Label(self, text="Input Format:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_format_var = tk.StringVar(self)
        self.input_format_var.set("fasta")
        self.input_format_menu = tk.OptionMenu(self, self.input_format_var, *self.supported_formats)
        self.input_format_menu.grid(row=1, column=1, padx=5, pady=5, sticky='w')

        # Create and place the convert button
        tk.Button(self, text="Open", command=self.open_file).grid(row=2, column=1, padx=5, pady=20)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    def disable_menu(self):
        pass

    def enable_menu(self):
        pass

    def browse_input_file(self):
        filename = filedialog.askopenfilename(title="Select input file")
        self.input_file_entry.delete(0, tk.END)
        self.input_file_entry.insert(0, filename)

    def open_file(self):
        self.in_file = self.input_file_entry.get()
        self.in_format = self.input_format_var.get()

        if not is_valid_path(self.in_file):
            TKMessageBox("Error", f"Please select a valid path for input file.", messagebox.showerror, self).show()
            return

        self.is_closed_succesfully = True
        self.destroy()


class SequenceSaveFileDialog(ModalWindow):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(self.parent, resizable=False)
        self.title("Sequence File Dialog")

        self.supported_formats = list(SeqIO._FormatToWriter.keys())

        #     "fasta": FastaIO.FastaWriter,
        #     "fasta-2line": FastaIO.FastaTwoLineWriter,
        #     "gb": InsdcIO.GenBankWriter,
        #     "genbank": InsdcIO.GenBankWriter,
        #     "embl": InsdcIO.EmblWriter,
        #     "imgt": InsdcIO.ImgtWriter,
        #     "nib": NibIO.NibWriter,
        #     "phd": PhdIO.PhdWriter,
        #     "pir": PirIO.PirWriter,
        #     "fastq": QualityIO.FastqPhredWriter,
        #     "fastq-sanger": QualityIO.FastqPhredWriter,
        #     "fastq-solexa": QualityIO.FastqSolexaWriter,
        #     "fastq-illumina": QualityIO.FastqIlluminaWriter,
        #     "qual": QualityIO.QualPhredWriter,
        #     "seqxml": SeqXmlIO.SeqXmlWriter,
        #     "sff": SffIO.SffWriter,
        #     "tab": TabIO.TabWriter,
        #     "xdna": XdnaIO.XdnaWriter,

        # init
        self.is_closed_succesfully = False
        self.out_file = None
        self.out_format = None
        self.selected_annots = None

        # Create and place labels and entry widgets
        tk.Label(self, text="Output File:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.output_file_entry = tk.Entry(self, width=30)
        self.output_file_entry.grid(row=0, column=1, padx=5, pady=5, sticky='we')

        tk.Button(self, text="Browse", command=self.browse_output_file).grid(row=0, column=2, padx=5, pady=5)

        tk.Label(self, text="Output Format:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.output_format_var = tk.StringVar(self)
        self.output_format_var.set(self.supported_formats[0])
        self.output_format_menu = tk.OptionMenu(self, self.output_format_var, *self.supported_formats)
        self.output_format_menu.grid(row=1, column=1, padx=5, pady=5, sticky='w')

        # Keep Annots Label Frame
        self.annots_frame = ttk.LabelFrame(self, text="Keep Annotations")
        self.annots_frame.grid(row=2, column=0, columnspan=4, padx=10, pady=5, sticky="ew")

        # Checkboxes
        self.name_var = tk.BooleanVar(value=True)
        self.descrp_var = tk.BooleanVar(value=True)
        self.annot_var = tk.BooleanVar(value=True)
        self.feature_var = tk.BooleanVar(value=True)
        self.letter_annot_var = tk.BooleanVar(value=True)
        self.dbxrefs_var = tk.BooleanVar(value=True)

        self.name_check = ttk.Checkbutton(self.annots_frame, text="Name", variable=self.name_var)
        self.descrp_check = ttk.Checkbutton(self.annots_frame, text="Description", variable=self.descrp_var)
        self.annot_check = ttk.Checkbutton(self.annots_frame, text="Annotation", variable=self.annot_var)
        self.feature_check = ttk.Checkbutton(self.annots_frame, text="Feature", variable=self.feature_var)
        self.letter_annot_check = ttk.Checkbutton(self.annots_frame, text="Letter Annotation",
                                                  variable=self.letter_annot_var)
        self.dbxrefs_check = ttk.Checkbutton(self.annots_frame, text="Database Cross-Reference",
                                             variable=self.dbxrefs_var)

        self.name_check.grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.descrp_check.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        self.annot_check.grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.feature_check.grid(row=1, column=1, sticky="w", padx=5, pady=2)
        self.letter_annot_check.grid(row=2, column=0, sticky="w", padx=5, pady=2)
        self.dbxrefs_check.grid(row=2, column=1, sticky="w", padx=5, pady=2)

        # Create and place the convert button
        tk.Button(self, text="Save", command=self.save_file).grid(row=3, column=1, padx=5, pady=20)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    def disable_menu(self):
        pass

    def enable_menu(self):
        pass

    def browse_output_file(self):
        filename = filedialog.asksaveasfilename(title="Select output file")
        self.output_file_entry.delete(0, tk.END)
        self.output_file_entry.insert(0, filename)

    def save_file(self):
        self.out_file = self.output_file_entry.get()
        directory = os.path.dirname(self.out_file)
        if not os.path.exists(directory):
            TKMessageBox("Error", f"Please select a valid path.", messagebox.showerror, self).show()
            return

        self.out_format = self.output_format_var.get()
        self.selected_annots = {
            'Name': self.name_var.get(),
            'Desc': self.descrp_var.get(),
            'Annot': self.annot_var.get(),
            'Feature': self.feature_var.get(),
            'Letter Annot': self.letter_annot_var.get(),
            'Dbxrefs': self.dbxrefs_var.get()
        }
        self.is_closed_succesfully = True
        self.destroy()

