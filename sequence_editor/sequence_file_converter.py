import os
import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import Phylo
from utils.model_window import ModalWindow
from Bio import AlignIO
from Bio import SeqIO

from utils.utils import is_valid_path
from utils.windows import TKMessageBox


class ConvertSequenceFormat(ModalWindow):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(self.parent, resizable=False)
        self.title("Sequence File Format Converter")

        self.in_file = None
        self.in_format = None
        self.out_file = None
        self.out_format = None
        self.is_closed_succesfully = False

        self.in_supported_formats = list(SeqIO._FormatToIterator.keys())

        list1 = list(SeqIO._FormatToString.keys())
        list2 = list(SeqIO._FormatToWriter.keys())
        self.out_supported_formats = []
        for x in list1 + list2:
            if x not in self.out_supported_formats:
                self.out_supported_formats.append(x)

        # Create and place labels and entry widgets
        tk.Label(self, text="Input File:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_file_entry = tk.Entry(self, width=30)
        self.input_file_entry.grid(row=0, column=1, padx=5, pady=5, sticky='we')

        tk.Button(self, text="Browse", command=self.browse_input_file).grid(row=0, column=2, padx=5, pady=5)

        tk.Label(self, text="Input Format:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.input_format_var = tk.StringVar(self)
        self.input_format_var.set(self.in_supported_formats[0])
        self.input_format_menu = tk.OptionMenu(self, self.input_format_var, *self.in_supported_formats)
        self.input_format_menu.grid(row=1, column=1, padx=5, pady=5, sticky='w')

        tk.Label(self, text="Output File:").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.output_file_entry = tk.Entry(self, width=30)
        self.output_file_entry.grid(row=2, column=1, padx=5, pady=5, sticky='we')

        tk.Button(self, text="Browse", command=self.browse_output_file).grid(row=2, column=2, padx=5, pady=5)

        tk.Label(self, text="Output Format:").grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
        self.output_format_var = tk.StringVar(self)
        self.output_format_var.set(self.out_supported_formats[0])
        self.output_format_menu = tk.OptionMenu(self, self.output_format_var, *self.out_supported_formats)
        self.output_format_menu.grid(row=3, column=1, padx=5, pady=5, sticky='w')

        # Create and place the convert button
        tk.Button(self, text="Convert", command=self.convert_file).grid(row=4, column=1, padx=5, pady=20)

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

    def browse_output_file(self):
        filename = filedialog.asksaveasfilename(title="Select output file")
        self.output_file_entry.delete(0, tk.END)
        self.output_file_entry.insert(0, filename)

    def convert_file(self):
        self.in_file = self.input_file_entry.get()
        self.in_format = self.input_format_var.get()

        if not is_valid_path(self.in_file):
            TKMessageBox("Error", f"Please select a valid path for input file.", messagebox.showerror, self).show()
            return

        self.out_file = self.output_file_entry.get()
        self.out_format = self.output_format_var.get()

        directory = os.path.dirname(self.out_file)
        if not os.path.exists(directory):
            TKMessageBox("Error", f"Please select a valid path for output file.", messagebox.showerror, self).show()
            return

        self.is_closed_succesfully = True
        self.destroy()


