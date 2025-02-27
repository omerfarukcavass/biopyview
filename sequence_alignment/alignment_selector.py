from utils.model_window import ModalWindow
import tkinter as tk
from tkinter import Spinbox


class SelectAlignment(ModalWindow):
    def __init__(self, parent, no_of_tables):
        self.parent = parent
        self.selected_table = 1  # starts from 1
        super().__init__(self.parent, resizable=False)

        # Create a new window for editing
        self.title("Select Alignment")

        # Create and place widgets
        self.label = tk.Label(self, text=f"This file consists of {no_of_tables} alignments. Select alignment to open:")
        self.label.pack(padx=10, pady=10)

        # self.spinbox = Spinbox(self, from_=1, to=no_of_tables)
        # self.spinbox.pack(padx=10, pady=10)

        self.selected_table_var = tk.StringVar(value = "1.alignment")
        self.aligns = [f"{i}.alignment" for i in range(1, no_of_tables + 1)]
        self.align_selector = tk.OptionMenu(self, self.selected_table_var, *self.aligns)
        self.align_selector.pack(padx=10, pady=10)

        # Add a "Save" button
        self.save_button = tk.Button(self, text="Confirm", command=self.save_changes)
        self.save_button.pack(padx=10, pady=10)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    # Function to save the new values and update the Treeview
    def save_changes(self):
        self.selected_table = int(self.selected_table_var.get().split(".")[0])
        self.destroy()
