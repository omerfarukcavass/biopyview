import os
import time
import tkinter as tk
from tkinter import ttk
from utils.model_window import ModalWindow
from utils.tooltip import Tooltip
from tkinter import Spinbox
from motif.motif import MotifSearch, FixedLengthMotif, VariableLengthMotif
from Bio.motifs import Motif
from utils.utils import is_dna, is_rna, is_amino_acid, DNA_BASES, RNA_BASES, AA_BASES, toggle_menu_on_off
from utils.utils import print_all_table_items
from tkinter import filedialog, messagebox
from utils.windows import MessageBox, TKMessageBox
from Bio import motifs
from viewer.sequence_viewer import SequenceViewer
from Bio.Seq import Seq
import threading

MOTIF_EMPTY_MESSAGE = "Motif does not have instances. Look at Counts table."


class MotifFinder(ModalWindow):
    def __init__(self, parent, motif_search: MotifSearch, seq_viewer:SequenceViewer):
        self.parent = parent  # main window, App class
        self.seq_viewer =  seq_viewer
        super().__init__(self.parent, grab_set=False, resizable=False)  # grab_set is False, user can select letters
        self.title("Motif Finder")

        # Init
        self.residue_opts = ["DNA", "RNA", "Aminoacid"]
        self.len_type_opts = ["Fixed Length", "Variable Length"]
        self.create_new_txt = "Create new"
        self.counts_table, self.pwm_table, self.pssm_table = None, None, None
        self.menu_bar = None
        self.menu_buttons = ["Edit Motif", "Edit Selected Row", "Complete Editing Motif", "Cancel Editing Motif", "Edit pseudo-counts", "Edit background", "Save motif"]
        self.motifs = None
        self.notebook_tabs = ["Instances", "Counts", "PWM", "PSSM", "Search Options"]
        self.motif_edit_on = False
        self.motif_edit_instance_mode = True
        self.search_manager = None
        # motifs and search options
        self.motif_search = motif_search
        self.current_motif_index = 0
        # save
        self.save_error = None
        self.message_box = None
        self.stop_search_lock = threading.Lock()

        # create GUI
        self.create_menu()
        self.create_widgets()
        self.show_curr_motif()

        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def on_close(self):
        # if search is not stopped
        if self.start_button["text"] != "Start":
            TKMessageBox("Error", "Motif search is on. Stop search before exiting. ", messagebox.showerror,
                         self).show()
            return

        quit = TKMessageBox("Warning", "Are you sure to quit now ?  Unsaved motifs will be lost.", messagebox.askyesno, self).show()
        if quit:
            self.parent.motif_finder_open = False
            self.destroy()
        else:
            return  # kapatmaktan vazgeçti

    def create_widgets(self):
        # Create the main frame
        self.main_frame = ttk.Frame(self)
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        self.create_motif_frame()
        self.create_center_frame()
        self.show_motif_and_center_frame(first_time=True)
        self.create_search_frame()

    def create_motif_frame(self):
        # Create a frame to settings
        self.motif_frame = tk.Frame(self.main_frame)
        self.motif_frame.pack(fill=tk.BOTH, expand=False)

        # Create a notebook in the settings frame
        motif_notebook = ttk.Notebook(self.motif_frame)
        motif_notebook.pack(fill=tk.BOTH, expand=True)

        # Create a tab for the settings buttons
        motif_tab = ttk.Frame(motif_notebook)
        motif_notebook.add(motif_tab, text='Motif')

        # Configure grid columns in the setting_frame
        motif_tab.columnconfigure(0, weight=1)
        motif_tab.columnconfigure(1, weight=1)
        motif_tab.columnconfigure(2, weight=1)
        motif_tab.columnconfigure(3, weight=1)
        motif_tab.rowconfigure(0, weight=1)

        # motif selector frame
        self.motif_selector_frame = ttk.Frame(motif_tab)
        self.motif_selector_frame.grid(row=0, column=0)

        # motif selector label
        motif_selector_label = ttk.Label(self.motif_selector_frame, text="Motif:")
        motif_selector_label.pack(side="top")

        # motif selector option menu
        self.selected_motif = tk.StringVar()
        size = max(1, len(self.motif_search.motifs))
        self.motifs = [f"{i}.motif" for i in range(1, size + 1)] + [self.create_new_txt]
        self.selected_motif.set(self.motifs[0])
        self.create_motif_selector()

        # residue type frame
        self.residue_type_frame = ttk.Frame(motif_tab)
        self.residue_type_frame.grid(row=0, column=1)

        # residue type label
        residue_type_label = ttk.Label(self.residue_type_frame, text="Residue:")
        residue_type_label.pack(side="top")

        # residue type option menu
        self.selected_residue = tk.StringVar()
        self.selected_residue.set(self.residue_opts[0])
        self.residue_selector = tk.OptionMenu(self.residue_type_frame, self.selected_residue, *self.residue_opts,
                                              command=self.toggle_residue)
        self.residue_selector.config(state=tk.DISABLED)
        self.residue_selector.pack(side="top")

        # length type frame
        self.type_frame = ttk.Frame(motif_tab)
        self.type_frame.grid(row=0, column=2)

        # length type label
        type_label = ttk.Label(self.type_frame, text="Length Type:")
        type_label.pack(side="top")

        # length type option menu
        self.selected_length_type = tk.StringVar()
        self.selected_length_type.set(self.len_type_opts[0])
        self.length_type_selector = tk.OptionMenu(self.type_frame, self.selected_length_type, *self.len_type_opts,
                                                  command=self.toggle_length_type)
        self.length_type_selector.config(state=tk.DISABLED)
        self.length_type_selector.pack(side="top")

        # length frame
        self.length_frame = ttk.Frame(motif_tab)
        self.length_frame.grid(row=0, column=3)

        # length label
        length_label = ttk.Label(self.length_frame, text="Length:")
        length_label.pack(side="top")

        # length spinbox
        self.selected_length = tk.IntVar()
        self.length_selector = Spinbox(self.length_frame, from_=1, to=10000, width=5, textvariable=self.selected_length,
                                       command=self.toggle_length)
        self.length_selector.config(state=tk.DISABLED)
        self.length_selector.pack(side="top")

    def create_center_frame(self):
        # Create a frame to hold notebook
        self.center_frame = tk.Frame(self.main_frame)

        # Create notebook
        self.notebook = ttk.Notebook(self.center_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Bind the <<NotebookTabChanged>> event
        self.notebook.bind('<<NotebookTabChanged>>', self.on_tab_change)

        # Add tabs to the Notebook
        self.instances_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.instances_tab, text=self.notebook_tabs[0])

        # instances text
        self.instances_text_widget = tk.Text(self.instances_tab, state="disabled", font=("Courier", 16), spacing1=2, spacing2=2,
                                             spacing3=2, height=10, wrap=tk.NONE)

        # Create a vertical scrollbar for the sequence names Text widget
        v_scrollbar = tk.Scrollbar(self.instances_tab, orient=tk.VERTICAL, command=self.instances_text_widget.yview)
        v_scrollbar.place(relx=1, rely=0, relheight=1, anchor='ne')
        self.instances_text_widget.configure(yscrollcommand=v_scrollbar.set)

        # Create a horizontal scrollbar for the Text widget
        h_scrollbar = tk.Scrollbar(self.instances_tab, orient=tk.HORIZONTAL,
                                                  command=self.instances_text_widget.xview)
        h_scrollbar.place(relx=0, rely=1, relwidth=1, anchor='sw')
        self.instances_text_widget.configure(xscrollcommand=h_scrollbar.set)

        self.instances_text_widget.pack(fill=tk.BOTH, expand=True)

        # Other tabs
        self.counts_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.counts_tab, text=self.notebook_tabs[1])

        self.pwm_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.pwm_tab, text=self.notebook_tabs[2])

        self.pssm_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.pssm_tab, text=self.notebook_tabs[3])

        # Create a tab for the search opts
        self.opt_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.opt_tab, text='Search Options')

        # Configure grid columns in the search_tab
        self.opt_tab.columnconfigure(0, weight=1)
        self.opt_tab.columnconfigure(1, weight=1)
        self.opt_tab.columnconfigure(2, weight=1)
        self.opt_tab.rowconfigure(0, weight=1)
        self.opt_tab.rowconfigure(1, weight=1)
        self.opt_tab.rowconfigure(2, weight=1)

        # Add search options to the search_tab
        options_frame = ttk.Frame(self.opt_tab)
        options_frame.grid(row=0, column=0, columnspan=3, pady=5, padx=5, sticky="nsew")
        options_frame.columnconfigure(0, weight=1)
        options_frame.columnconfigure(1, weight=1)
        options_frame.rowconfigure(0, weight=1)

        self.search_scope_var = tk.StringVar(value="all")
        self.search_start_var = tk.StringVar(value="beginning")
        self.search_type_var = tk.StringVar(value="exact")

        # Search Scope
        scope_frame = ttk.LabelFrame(options_frame, text="Search Scope")
        scope_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        # all
        self.search_all_radio = ttk.Radiobutton(
            scope_frame, text="Search in all sequences", variable=self.search_scope_var, value="all"
        )
        self.search_all_radio.pack(anchor="w", padx=5, pady=2)
        self.search_all_radio.config(command=self.toggle_start_options)

        # selected
        self.search_selected_radio = ttk.Radiobutton(
            scope_frame, text="Search in selected sequences", variable=self.search_scope_var, value="selected"
        )
        self.search_selected_radio.pack(anchor="w", padx=5, pady=2)
        self.search_selected_radio.config(command=self.toggle_start_options)

        # Search Starting Point
        thresold_frame = ttk.LabelFrame(options_frame, text="Search Starting Point")
        thresold_frame.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

        # beginning
        self.start_beginning_radio = ttk.Radiobutton(
            thresold_frame, text="Start from beginning", variable=self.search_start_var, value="beginning"
        )
        self.start_beginning_radio.pack(anchor="w", padx=5, pady=2)
        # self.start_beginning_radio.config(state=tk.DISABLED)

        # selected
        self.start_selected_radio = ttk.Radiobutton(
            thresold_frame, text="Start from selected sequence(s)", variable=self.search_start_var, value="selected"
        )
        self.start_selected_radio.pack(anchor="w", padx=5, pady=2)
        # self.start_selected_radio.config(state=tk.DISABLED)

        # Search Type
        type_frame = ttk.LabelFrame(options_frame, text="Search Type")
        type_frame.grid(row=1, column=0, columnspan=1, padx=5, pady=5, sticky="nsew")

        self.exact_search_radio = ttk.Radiobutton(
            type_frame, text="Exact search", variable=self.search_type_var, value="exact", command=self.toggle_threshold
        )
        self.exact_search_radio.pack(anchor="w", padx=5, pady=2)

        self.pssm_search_radio = ttk.Radiobutton(
            type_frame, text="PSSM search", variable=self.search_type_var, value="pssm", command=self.toggle_threshold
        )
        self.pssm_search_radio.pack(anchor="w", padx=5, pady=2)

        # PSSM threshold
        thresold_frame = ttk.LabelFrame(options_frame, text="PSSM Threshold")
        thresold_frame.grid(row=1, column=1, padx=5, pady=5, sticky="nsew")

        self.threshold_entry = ttk.Entry(thresold_frame, width=5)
        self.threshold_entry.pack(anchor="w", padx=5, pady=2)
        self.threshold_entry.config(state=tk.DISABLED)

    def create_motif_selector(self):
        self.motif_selector = tk.OptionMenu(self.motif_selector_frame, self.selected_motif, *self.motifs,
                                            command=self.motif_selected)
        self.motif_selector.pack(side="top")

    def create_search_frame(self):
        # Create the search frame
        self.search_frame = tk.Frame(self.main_frame)
        self.search_frame.pack(fill=tk.BOTH, expand=False)

        # Create a notebook in the search frame
        search_notebook = ttk.Notebook(self.search_frame)
        search_notebook.pack(fill=tk.BOTH, expand=True)

        # Create a tab for the search buttons
        search_tab = ttk.Frame(search_notebook)
        search_notebook.add(search_tab, text='Search')

        # Configure grid columns in the search_tab
        search_tab.columnconfigure(0, weight=1)
        search_tab.columnconfigure(1, weight=1)
        search_tab.columnconfigure(2, weight=1)
        search_tab.columnconfigure(3, weight=1)
        search_tab.rowconfigure(0, weight=1)
        search_tab.rowconfigure(1, weight=1)

        # Add buttons to the search_tab
        self.start_button = ttk.Button(search_tab, text="Start", command=self.toggle_start_button)
        self.start_button.grid(row=0, column=0, padx=5, pady=5)

        self.prev_button = ttk.Button(search_tab, text="◀", width=2, command=self.prev_match)
        self.prev_button.grid(row=0, column=1, padx=5, pady=5)
        self.prev_button.config(state=tk.DISABLED)

        self.next_button = ttk.Button(search_tab, text="▶", width=2, command=self.next_match)
        self.next_button.grid(row=0, column=2, padx=5, pady=5)
        self.next_button.config(state=tk.DISABLED)

        self.min_max_button = ttk.Button(search_tab, text="Minimize", command=self.toggle_screen)
        self.min_max_button.grid(row=0, column=3, padx=5, pady=5)
        # tt = Tooltip(self.min_max_button, "Minimize screen.")

        self.status_bar = tk.Label(search_tab, text="-")
        self.status_bar.grid(row=1, column=0, columnspan=4, padx=5, pady=5, sticky=tk.W)

    def create_menu(self):
        # Create the menu
        self.menu_bar = tk.Menu(self)
        self.config(menu=self.menu_bar)

        self.file_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.menu_bar.add_cascade(label="File", menu=self.file_menu)

        self.file_menu.add_command(label=self.menu_buttons[0], command=self.menu_edit_motif)
        self.file_menu.add_command(label=self.menu_buttons[1], command=self.menu_edit_selected_row)
        self.file_menu.add_command(label=self.menu_buttons[2], command=self.menu_complete_editing)
        self.file_menu.add_command(label=self.menu_buttons[3], command=self.menu_cancel_editing)
        self.file_menu.add_command(label=self.menu_buttons[4], command=self.menu_edit_pseudo_counts)
        self.file_menu.add_command(label=self.menu_buttons[5], command=self.menu_edit_background)
        self.file_menu.add_command(label=self.menu_buttons[6], command=self.menu_save_motif)

        self.file_menu.entryconfig(1, state=tk.DISABLED)  # "Edit Selected Row"
        self.file_menu.entryconfig(2, state=tk.DISABLED)  # "Complete Editing Motif"
        self.file_menu.entryconfig(3, state=tk.DISABLED)  # "Cancel Editing Motif"

    def show_curr_motif(self, update_motif_settings=True, update_instances_tab=True, update_counts_tab=True, update_pwm_tab=True, update_pssm_tab=True, update_search_tab=True):
        """
        Generic method for displaying a motif in screen.

        :param update_motif_settings:
        :param update_instances_tab:
        :param update_counts_tab:
        :param update_pwm_tab:
        :param update_pssm_tab:
        :return:
        """
        if len(self.motif_search.motifs) == 0:
            return

        if len(self.motif_search.motifs) == 1 and self.motif_search.motifs[0] is None:
            return

        curr_motif = self.motif_search.motifs[self.current_motif_index]

        # If not a fixed motif, disable pseudo and bg menu
        if not isinstance(curr_motif, Motif):
            self.file_menu.entryconfig(4, state=tk.DISABLED)  # pseudo
            self.file_menu.entryconfig(5, state=tk.DISABLED)  # bg
        else:
            self.file_menu.entryconfig(4, state=tk.NORMAL)  # pseudo
            self.file_menu.entryconfig(5, state=tk.NORMAL)  # bg

        if update_motif_settings:
            # set residue type
            residue_type = self.residue_opts[0] if curr_motif.alphabet == DNA_BASES else (
                self.residue_opts[1] if curr_motif.alphabet == RNA_BASES else self.residue_opts[2])
            self.selected_residue.set(residue_type)

            # set length type
            is_fixed = isinstance(curr_motif, Motif)
            self.selected_length_type.set(self.len_type_opts[0] if is_fixed else self.len_type_opts[1])

            # set length
            if is_fixed:
                self.selected_length.set(len(curr_motif))
            else:
                self.toggle_length_type()  # disable spinbox

        if update_instances_tab:
            # update text widget
            if curr_motif.instances is None:
                self.instances_text_widget.configure(state="normal")
                self.instances_text_widget.delete("1.0", tk.END)
                self.instances_text_widget.update_idletasks()
                self.instances_text_widget.insert(tk.END, MOTIF_EMPTY_MESSAGE + "\n")
                self.instances_text_widget.configure(state="disabled")
            else:
                self.instances_text_widget.configure(state="normal")
                self.instances_text_widget.delete("1.0", tk.END)
                self.instances_text_widget.update_idletasks()

                for seq in curr_motif.instances:
                    self.instances_text_widget.insert(tk.END, seq + "\n")

                self.instances_text_widget.configure(state="disabled")

        if update_counts_tab:
            # delete pfm, pwm, pssm tabs
            if self.counts_table is not None:
                self.counts_table.destroy()
            self.counts_table = None

            # update tab if it is a "fixed" motif.
            if isinstance(curr_motif, Motif):

                # create pfm
                self.counts_table = self.create_treeview(self.counts_tab, ("Position", *curr_motif.alphabet))
                values_generator = ((i+1, *pos_counts) for i, (*pos_counts,) in enumerate(zip(*curr_motif.counts.values())))
                self.populate_treeview(self.counts_table, values_generator)
                self.counts_table.pack(expand=True, fill="both")
                self.notebook.tab(1, state=tk.NORMAL)
            else:
                self.notebook.tab(1, state=tk.DISABLED)

        if update_pwm_tab:
            # delete pfm, pwm, pssm tabs
            if self.pwm_table is not None:
                self.pwm_table.destroy()
            self.pwm_table = None

            # update tab if it is a "fixed" motif.
            if isinstance(curr_motif, Motif):

                # create pwm
                self.pwm_table = self.create_treeview(self.pwm_tab, ("Position", *curr_motif.alphabet))
                values_generator = ((i + 1, *pos_counts) for i, (*pos_counts,) in
                                    enumerate(zip(*curr_motif.pwm.values())))
                self.populate_treeview(self.pwm_table, values_generator)
                self.pwm_table.pack(expand=True, fill="both")
                self.notebook.tab(2, state=tk.NORMAL)

            else:
                self.notebook.tab(2, state=tk.DISABLED)

        if update_pssm_tab:
            # delete pfm, pwm, pssm tabs
            if self.pssm_table is not None:
                self.pssm_table.destroy()
            self.pssm_table = None

            # update tab if it is a "fixed" motif.
            if isinstance(curr_motif, Motif):

                # create pssm
                self.pssm_table = self.create_treeview(self.pssm_tab, ("Position", *curr_motif.alphabet))
                values_generator = ((i + 1, *pos_counts) for i, (*pos_counts,) in
                                    enumerate(zip(*curr_motif.pssm.values())))
                self.populate_treeview(self.pssm_table, values_generator)
                self.pssm_table.pack(expand=True, fill="both")
                self.notebook.tab(3, state=tk.NORMAL)

            else:
                self.notebook.tab(3, state=tk.DISABLED)

        if update_search_tab:
            if not isinstance(curr_motif, Motif):  # if var len, pssm disabled.
                self.pssm_search_radio.config(state=tk.DISABLED)
            else:
                self.pssm_search_radio.config(state=tk.NORMAL)

    def create_treeview(self, tab, column_names):
        colns = [str(i) for i in range(len(column_names))]
        tree = ttk.Treeview(tab, show='headings', selectmode='browse', columns=colns)

        for cid, column_name in zip(colns, column_names):
            tree.heading(cid, text=column_name)
            tree.column(cid, width=50)

        # Create vertical scrollbar
        vscrollbar = ttk.Scrollbar(tab, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=vscrollbar.set)
        vscrollbar.place(relx=1, rely=0, relheight=1, anchor='ne')

        # Create horizontal scrollbar
        hscrollbar = ttk.Scrollbar(tab, orient="horizontal", command=tree.xview)
        tree.configure(xscrollcommand=hscrollbar.set)
        hscrollbar.place(relx=0, rely=1, relwidth=1, anchor='sw')

        # Bind the Configure event to adjust column widths
        #tree.bind("<Configure>", lambda event: self.tree_on_configure(event, tree))

        return tree

    def populate_treeview(self, tree, data_iterator):
        for values in data_iterator:
            tree.insert("", index="end", values=values)

    def clear_treeview(self, tree, last_n_rows=None):
        if last_n_rows is None:  # clear all rows
            for item in tree.get_children():
                tree.delete(item)
        else:  # clear last n rows
            for item in tree.get_children()[-last_n_rows:]:
                tree.delete(item)

    def extend_treeview(self, tree, n_row):
        row_size = len(tree['columns']) - 1
        last_index = int(tree.item(tree.get_children()[-1], "values")[0]) + 1 if len(tree.get_children()) > 0 else 1

        for i in range(last_index, last_index + n_row):
            tree.insert("", index="end", values=[i] + [0]*row_size)

    def rename_treeview_columns(self, tree, new_columns):
        cids = list(tree['columns'])[1:]  # cid for letter columns
        if len(cids) != len(new_columns):
            print("Column number does not match !")
            return

        for cid, new_name in zip(cids, new_columns):
            tree.heading(cid, text=new_name)

    def tree_on_configure(self, event, tree):
        # Adjust column widths when the window size changes
        for column in tree['columns']:
            tree.column(column, width=50, minwidth=50)

    def toggle_start_button(self):
        if self.start_button["text"] == "Start":
            self.start_searching()
        else:
            self.stop_searching()

    def start_searching(self):
        # get opts
        search_in_all = self.search_scope_var.get() == "all"
        start_from_selected = search_in_all and self.search_start_var.get() == "selected"
        pssm_search = self.search_type_var.get() == "pssm"
        pssm_threshold_s = self.threshold_entry.get()

        if not search_in_all and not self.seq_viewer.any_selection(self.seq_viewer.seq_text_widget):
            TKMessageBox("Error", f"There is no selected sequence to search!", messagebox.showerror, self).show()
            return

        if start_from_selected and not self.seq_viewer.any_selection(self.seq_viewer.seq_text_widget):
            TKMessageBox("Error", f"There is no selected sequence to start from!", messagebox.showerror, self).show()
            return

        curr_search_motif = self.motif_search.motifs[self.current_motif_index]
        if curr_search_motif.instances is None and not pssm_search:
            TKMessageBox(
                "Error",
                f"This motif does not have instances, so Exact Search cannot be done! \nYou might switch to PSSM search.",
                messagebox.showerror,
                self
            ).show()
            return

        if curr_search_motif.alphabet is not DNA_BASES and pssm_search:
            TKMessageBox("Error", f"PSSM search is only available for DNA motifs!", messagebox.showerror, self).show()
            return

        # threshold float olmazsa meesagebox
        if pssm_search:
            try:
                pssm_threshold = float(pssm_threshold_s)
            except ValueError:
                TKMessageBox("Error", f"PSSM threshold must be a float.", messagebox.showerror, self).show()
                return
        else:
            pssm_threshold = None

        if pssm_search and pssm_threshold <0:
            TKMessageBox("Error", f"PSSM threshold must be a non-negative value.", messagebox.showerror, self).show()
            return

        # disable other GUI elements
        self.update_search_gui(True)
        self.start_button.config(text="Stop", style="Accent.TButton")

        # define search area
        if not search_in_all:  # search in selected
            pos_1, pos_2, pos_3, pos_4 = self.seq_viewer.get_selection_bounds()
            search_start_line_index = int(pos_1.split('.')[0])
            search_stop_line_index = int(pos_3.split('.')[0])

        elif start_from_selected:
            pos_1, pos_2, pos_3, pos_4 = self.seq_viewer.get_selection_bounds()
            search_start_line_index = int(pos_1.split('.')[0])
            search_stop_line_index = self.seq_viewer.get_no_of_lines(self.seq_viewer.seq_text_widget)-1

        else:
            search_start_line_index = self.seq_viewer.start_offset_in_text_widget
            search_stop_line_index = self.seq_viewer.get_no_of_lines(self.seq_viewer.seq_text_widget)-1

        # create search data (we don't call update_data in seq_viewer as it cost more than this. )
        self.status_bar.config(text="Preparing...")
        self.status_bar.update_idletasks()
        search_data = []

        for line_index in range(search_start_line_index, search_stop_line_index + 1):
            seq_name_line = self.seq_viewer.get_line(self.seq_viewer.seq_names_text_widget, line_index)
            seq_id = self.seq_viewer.parse_seq_id(seq_name_line)
            if seq_id in self.seq_viewer.edited_seq_ids:
                seq_line = self.seq_viewer.get_line(self.seq_viewer.seq_text_widget, line_index)
            else:
                seq_line = self.seq_viewer.data.get_record(seq_id).seq
            search_data.append(seq_line)

        if len(search_data) == 0:
            TKMessageBox("Error", "No sequence available for search !", messagebox.showinfo, self).show()
            self.stop_searching_gui()
            return

        self.search_manager = SearchManager(self.seq_viewer, self, curr_search_motif, pssm_search, pssm_threshold, search_data, search_start_line_index, search_stop_line_index)
        self.next_match()

    def prev_match(self):
        match = self.search_manager.prev()
        self.update_status_bar()
        self.seq_viewer.see_position(line_number=match.line_index, col_number=match.pos, highlight_len=match.length)
        self.seq_viewer.display_index_line()
        self.seq_viewer.paint_visible_lines_all()
        self.update_nav_buttons_status()

    def next_match(self):
        # disable prev and next if searching continues
        if not self.search_manager.search_completed:
            self.status_bar.config(text="Searching...")
            self.status_bar.update_idletasks()
            self.disable_nav_buttons()

        match = self.search_manager.next()  # search next item in a thread

        if match is not None:  # match already found, no query thread created
            self.show_current_match(match)

    def process_query_result(self):
        """
        This is called by query thread after query finishes.
        :return:
        """

        # --- Edge case: -----------------------
        # User clicked "stop" button but query thread had already got search_is_on as True, so it continued
        # and found a new match without realizing that user stopped search. Event order:
        # 1. query thread got search_is_on as True
        # 2. CPU switches to main thread and user click STOP button. Then CPU switch back to query thread.
        # 3. query thread found a match, thinking search is going on. It called this fun with found match.

        with self.stop_search_lock:
            if not self.search_manager.search_is_on:
                self.stop_searching_gui()
        # --------------------------------------

        # check if query thread find a new match or not.
        if self.search_manager.query_result is not None:  # new match found
            self.show_current_match(self.search_manager.query_result)
            self.search_manager.query_result = None  # reset for new query
        else:  # no match left
            self.update_status_bar()
            no_of_match = len(self.search_manager.search_matches)
            TKMessageBox("Info", "Motif search has completed..", messagebox.showinfo, self).show()

            if no_of_match == 0:
                self.stop_searching_gui()
                return

            self.update_nav_buttons_status()

    def show_current_match(self, match):
        self.update_status_bar()
        self.seq_viewer.see_position(line_number=match.line_index, col_number=match.pos, highlight_len=match.length)
        self.seq_viewer.display_index_line()
        self.seq_viewer.paint_visible_lines_all()
        self.update_nav_buttons_status()

    def update_status_bar(self):
        match_count = self.search_manager.get_match_count()
        if match_count is not None:
            self.status_bar.config(text=f"{self.search_manager.curr_match_index+1} of all matches.({match_count})")
        else:
            self.status_bar.config(text=f"{self.search_manager.curr_match_index+1} of all matches.")

    def update_nav_buttons_status(self):
        print("curr_match_index", self.search_manager.curr_match_index)
        # prev
        if self.search_manager.curr_match_index >= 1:
            self.prev_button.config(state=tk.NORMAL)
        else:
            self.prev_button.config(state=tk.DISABLED)

        # next
        if self.search_manager.search_completed:
            if self.search_manager.curr_match_index < len(self.search_manager.search_matches)-1:
                self.next_button.config(state=tk.NORMAL)
            else:
                self.next_button.config(state=tk.DISABLED)
        else:
            self.next_button.config(state=tk.NORMAL)

    def disable_nav_buttons(self):
        self.prev_button.config(state=tk.DISABLED)
        self.next_button.config(state=tk.DISABLED)

    def stop_searching(self):
        # if it is searching now, stops the search thread with search_is_on flag
        if self.search_manager.query_thread.is_alive():
            print("alive")
            self.status_bar.config(text="Stopping...")
            with self.stop_search_lock:
                self.search_manager.search_is_on = False
        else:
            print("dead")
            self.stop_searching_gui()

    def stop_searching_gui(self):
        """
        Run by main thread if query thread is NOT running now.
        Run by query thread after finishes its job if it is running now.
        :return:
        """
        # GUI elements update after stop
        self.status_bar.config(text="-")
        self.update_search_gui(False)
        self.start_button.config(text="Start", style="TButton")

    def update_search_gui(self, start_search: bool):
        first_state = tk.DISABLED if start_search else tk.NORMAL
        second_state = tk.DISABLED if not start_search else tk.NORMAL

        self.motif_selector.config(state=first_state)
        self.search_selected_radio.config(state=first_state)
        self.search_all_radio.config(state=first_state)
        self.exact_search_radio.config(state=first_state)
        self.pssm_search_radio.config(state=first_state)

        if self.search_scope_var.get() == "all":
            self.start_selected_radio.config(state=first_state)
            self.start_beginning_radio.config(state=first_state)

        if self.search_type_var.get() == "pssm":
            self.threshold_entry.config(state=first_state)

        self.file_menu.entryconfig(0, state=first_state)  # edit motif
        self.file_menu.entryconfig(4, state=first_state)  # edit pseudo
        self.file_menu.entryconfig(5, state=first_state)  # edit bg
        self.file_menu.entryconfig(6, state=first_state)  # save

        if start_search:  # prev off, next on
            self.prev_button.config(state=first_state)
            self.next_button.config(state=second_state)
        else:  # prev off, next off
            self.prev_button.config(state=second_state)
            self.next_button.config(state=second_state)

    def toggle_screen(self):
        if self.min_max_button["text"] == "Minimize":
            self.min_max_button.config(text="Maximize", style="Accent.TButton")
            self.hide_setting_and_center_frame()
        else:
            self.min_max_button.config(text="Minimize", style="TButton")
            self.show_motif_and_center_frame()

    # Function to hide
    def hide_setting_and_center_frame(self):
        self.motif_frame.pack_forget()
        self.center_frame.pack_forget()

    # Function to show
    def show_motif_and_center_frame(self, first_time=False):
        if first_time:
            self.motif_frame.pack(fill=tk.BOTH, expand=True)
            self.center_frame.pack(fill=tk.BOTH, expand=True)
        else:
            self.motif_frame.pack(fill=tk.BOTH, expand=True, before = self.search_frame)
            self.center_frame.pack(fill=tk.BOTH, expand=True, before = self.search_frame)

    def motif_selected(self, *args):
        """
        Switch btw motifs.
        :param args:
        :return:
        """
        motif_txt = self.selected_motif.get()

        # Create new motif
        if motif_txt == self.create_new_txt:

            # ask user edit mode: instances or counts
            w = SelectEditMode(self)

            if w.is_closed_succesfully:
                self.motif_edit_instance_mode = w.is_instance_mode

                # take last motif index, increment 1
                last_motif = int(self.motifs[-2].split(".")[0])
                new_motif= f"{last_motif+1}.motif"
                self.motifs.insert(-1, new_motif)

                # create "None" motif and append
                self.current_motif_index = last_motif
                self.motif_search.motifs.append(None)  # none motif

                # Delete OptionMenu and create new one
                self.motif_selector.destroy()
                self.create_motif_selector()
                self.selected_motif.set(new_motif)

                # Reset motif settings buttons
                self.selected_residue.set(self.residue_opts[0])
                self.selected_length_type.set(self.len_type_opts[0])
                self.selected_length.set(1)

                # reset text widget
                self.instances_text_widget.configure(state=tk.NORMAL)
                self.instances_text_widget.delete("1.0", tk.END)
                self.instances_text_widget.update_idletasks()
                self.instances_text_widget.configure(state=tk.DISABLED)

                # reset counts table
                self.recreate_counts_table(from_scratch=True)

                self.start_editing()

        else:
            # update current motif index to selected one
            new_index = int(motif_txt.split(".")[0])-1
            self.current_motif_index = new_index
            self.show_curr_motif()

    def menu_edit_motif(self, is_instance_mode=None):
        """
        Make current motif editable.
        :return:
        """

        # ask user edit mode: instances or counts
        w = SelectEditMode(self)

        if w.is_closed_succesfully:
            self.motif_edit_instance_mode = w.is_instance_mode
            self.start_editing()

    def start_editing(self):
        # make selected tab editable
        self.motif_edit_on = True
        self.notebook.select(self.instances_tab if self.motif_edit_instance_mode else self.counts_tab)
        self.switch_motif_editable(True)

    def create_motif(self, is_instance_mode):
        self.motif_edit_instance_mode = is_instance_mode
        self.recreate_counts_table()  # reset counts table

        # make selected tab editable
        self.motif_edit_on = True
        self.notebook.select(self.instances_tab if self.motif_edit_instance_mode else self.counts_tab)
        self.switch_motif_editable(True)

        self.file_menu.entryconfig(3, state=tk.DISABLED)  # "Cancel Editing Motif"

    def menu_edit_selected_row(self):
        EditSelectedRow(self, self.counts_table, ["Position"])

    def menu_complete_editing(self):

        # get residue type
        is_dna_selected = self.selected_residue.get() == self.residue_opts[0]
        is_rna_selected = self.selected_residue.get() == self.residue_opts[1]
        is_aa_selected = self.selected_residue.get() == self.residue_opts[2]
        alphabet = DNA_BASES if is_dna_selected else (RNA_BASES if is_rna_selected else AA_BASES)

        if self.motif_edit_instance_mode:
            # check if at least one letter exist
            seqs = self.get_instances_text_widget(to_upper=True)
            if len(seqs) == 0:
                print("No seqs exist.")
                TKMessageBox("Error", f"No sequence exist !", messagebox.showerror,
                             self).show()
                return

            # check if settings and instances are consistent
            # check residue type
            if is_dna_selected and any(not is_dna(seq) for seq in seqs):
                print("Non-DNA letter exist.")
                TKMessageBox("Error", f"Non-DNA letter exists!", messagebox.showerror, self).show()
                return

            if is_rna_selected and any(not is_rna(seq) for seq in seqs):
                print("Non-RNA letter exist.")
                TKMessageBox("Error", f"Non-RNA letter exists!", messagebox.showerror, self).show()
                return

            if is_aa_selected and any(not is_amino_acid(seq) for seq in seqs):
                print("Non-aminoacid letter exist.")
                TKMessageBox("Error", f"Non-aminoacid letter exists!", messagebox.showerror, self).show()
                return

            # check length type and size
            is_fixed_selected = self.selected_length_type.get() == self.len_type_opts[0]
            size = self.selected_length.get()
            if is_fixed_selected and any(len(seq)!=size for seq in seqs):
                print("Motif instances are not in selected length.")
                TKMessageBox("Error", f"Motif instances are not in selected length!", messagebox.showerror, self).show()
                return

            # All checks are passed. Create motif object
            if is_fixed_selected:
                motif = FixedLengthMotif(alphabet=alphabet, instances=seqs)
            else:
                motif = VariableLengthMotif(alphabet=alphabet, instances=seqs)

        else:
            # check if row sums are equal and not zero
            # first row
            all_items = self.counts_table.get_children()
            first_row_values = self.counts_table.item(all_items[0], "values")
            first_row_sum = sum(map(float, first_row_values[1:]))
            if first_row_sum == 0:
                print("Row sum cannot be zero !")
                TKMessageBox("Error", f"Row sum cannot be zero!", messagebox.showerror, self).show()
                return
            # other rows
            for row_id in all_items[1:]:
                row_values = self.counts_table.item(row_id, "values")
                row_sum = sum(map(float, row_values[1:]))
                if row_sum != first_row_sum:
                    print("Row sums are not equal !", row_sum, first_row_sum)
                    TKMessageBox("Error", f"Row sums are not equal!", messagebox.showerror, self).show()
                    return
                if row_sum == 0:
                    print("Row sum cannot be zero !")
                    TKMessageBox("Error", f"Row sum cannot be zero!", messagebox.showerror, self).show()
                    return

            # All checks are passed. Create motif object and update list.
            counts = self.create_counts_from_table(self.counts_table)
            motif = FixedLengthMotif(alphabet=alphabet, counts=counts)

        self.motif_search.motifs[self.current_motif_index] = motif

        # update motif tab editable
        self.switch_motif_editable(False)
        self.motif_edit_on = False

        # update notebook tabs
        if self.motif_edit_instance_mode:
            self.show_curr_motif(update_motif_settings=False, update_instances_tab=False)
        else:  # counts mode
            self.show_curr_motif(update_motif_settings=False, update_counts_tab=False)

    def menu_cancel_editing(self):
        # update motif tab editable
        self.switch_motif_editable(False)
        self.motif_edit_on = False

        # if "create new" clicked in motif selector and then this is called
        curr_motif = self.motif_search.motifs[self.current_motif_index]
        if curr_motif is None:
            # take last motif index, increment 1
            self.motifs.pop(-2)

            # remove "None" motif from list
            self.current_motif_index -= 1
            self.motif_search.motifs.pop(-1)

            # Delete OptionMenu and create new one
            self.motif_selector.destroy()
            self.create_motif_selector()
            self.selected_motif.set(self.motifs[-2])

        # update notebook tabs
        self.show_curr_motif()

    def menu_edit_pseudo_counts(self):
        # get new PseudoCounts and update motif
        motif = self.motif_search.motifs[self.current_motif_index]
        EditPseudoCounts(self, motif)

        # update pwm and pssm tabs
        self.show_curr_motif(update_motif_settings=False, update_instances_tab=False, update_counts_tab=False)

    def menu_edit_background(self):
        # get new background and update motif
        motif = self.motif_search.motifs[self.current_motif_index]
        EditBackground(self, motif)

        # update pwm and pssm tabs
        self.show_curr_motif(update_motif_settings=False, update_instances_tab=False, update_counts_tab=False, update_pwm_tab=False)

    def menu_save_motif(self):
        # filter fixed length motifs
        fixed_motifs = {f"{i+1}.motif":motif for i, motif in enumerate(self.motif_search.motifs) if isinstance(motif, Motif)}

        if len(fixed_motifs) > 0:
            smd = SaveMotifDialog(self, motifs=fixed_motifs)
            if smd.is_closed_succesfully:
                selected_motifs = smd.selected_motifs
                selected_format = smd.selected_format
                selected_path = smd.selected_path
                self.save_error = None

                # Create and start thread for writing
                thread = threading.Thread(target=self.save_motif_file_t, args=(selected_motifs, selected_format, selected_path))
                thread.start()

                # Popup window
                self.message_box = MessageBox(self, "Saving file...", window_title="Info", block_parent=False,
                                              ok_button=False)

                # Check write thread periodically
                self.check_save_thread(thread)
        else:
            TKMessageBox("Error", "No \"Fixed Length\" motif found to save.", messagebox.showerror, self).show()

    def save_motif_file_t(self, selected_motifs, selected_format, selected_path):  # thread fun
        # time.sleep(5)
        # Write data to file with selected format.
        try:
            with open(selected_path, "w") as out_handle:
                text = motifs.write(selected_motifs, selected_format)
                out_handle.write(text)
                print("Motifs written to selected_path.")
        except Exception as e:
            print(e)
            self.save_error = f"Error: {e}"
        else:
            self.save_error = None

    def check_save_thread(self, thread):
        # print("check_save_thread")
        if not thread.is_alive():
            self.message_box.destroy()
            if self.save_error is not None:  # thread is dead, so no lock needed
                TKMessageBox("Error", f"File couldn't be saved ! \n {self.save_error}", messagebox.showinfo, self).show()
            else:
                TKMessageBox("Info", "File saved successfully.", messagebox.showinfo, self).show()
        else:
            # Schedule the next check
            self.after(100, lambda: self.check_save_thread(thread))

    def create_counts_from_table(self, treeview):
        # Get letters
        cids = treeview["columns"][1:]
        colns = [treeview.heading(cid)["text"] for cid in cids]

        # Initialize a dictionary to store column values
        counts = {col: [] for col in colns}

        # Get all row identifiers
        rows = treeview.get_children()

        # Iterate through each row and gather values
        for row_id in rows:
            row_values = treeview.item(row_id, "values")
            for col, value in zip(colns, row_values[1:]):
                counts[col].append(float(value))

        return counts

    def get_instances_text_widget(self, to_upper=False):
        text_content = self.instances_text_widget.get("1.0", tk.END)
        if to_upper:
            seqs = [line.upper() for line in text_content.splitlines() if line.strip()]
        else:
            seqs = [line for line in text_content.splitlines() if line.strip()]
        return seqs

    def switch_motif_editable(self, enable):
        """
        Make motif editable or not based on edit mode (Instance or Counts).

        :param enable:
        :return:
        """
        first_state = tk.NORMAL if enable else tk.DISABLED  # selected state
        other_state = tk.NORMAL if not enable else tk.DISABLED  # opposite

        # motif selector
        self.motif_selector.config(state=other_state)

        # motif settings
        self.residue_selector.config(state=first_state)
        if self.motif_edit_instance_mode:  # Don't change length type selector if Counts mode
            self.length_type_selector.config(state=first_state)
        if first_state == tk.DISABLED or (first_state == tk.NORMAL and self.selected_length_type.get() == self.len_type_opts[0]):
            self.length_selector.config(state=first_state)

        # Instances or Counts
        if self.motif_edit_instance_mode:
            #  text widget
            self.instances_text_widget.configure(state=first_state)

            # delete MOTIF_EMPTY_MESSAGE if editing enabled
            if first_state == tk.NORMAL:
                instances = self.get_instances_text_widget()
                if len(instances) == 1 and instances[0] == MOTIF_EMPTY_MESSAGE:
                    self.instances_text_widget.delete("1.0", tk.END)
                    self.instances_text_widget.update_idletasks()
        else:  # counts mode
            # select_first_row
            self.select_first_row(self.counts_table)

            # adjust menu button
            self.file_menu.entryconfig(1, state=first_state)

        # search buttons
        self.start_button.config(state=other_state)
        # self.prev_button.config(state=other_state)
        # self.next_button.config(state=other_state)
        self.min_max_button.config(state=other_state)

        # tab switch enable/disable
        self.notebook.tab(1 if self.motif_edit_instance_mode else 0, state=other_state)
        self.notebook.tab(2, state=other_state)  # pwm
        self.notebook.tab(3, state=other_state)  # pssm
        self.notebook.tab(4, state=other_state)  # search opts

        #  menu
        self.file_menu.entryconfig(0, state=other_state)  # edit motif
        self.file_menu.entryconfig(2, state=first_state)  # complete editing
        self.file_menu.entryconfig(3, state=first_state)  # cancel editing
        self.file_menu.entryconfig(4, state=other_state)  # pseudo
        self.file_menu.entryconfig(5, state=other_state)  # bg
        self.file_menu.entryconfig(6, state=other_state)  # save


    def select_first_row(self, table):
        first_row_id = table.get_children()[0]
        table.selection_set(first_row_id)
        table.focus(first_row_id)

    def change_menu_edit_button(self, new_label):
        """
        Change menu button label
        :param index:
        :return:
        """

        # Delete label in menu
        if new_label is None:
            self.file_menu.delete(0)
        # Change label in menu
        else:
            # Get the index of the last command
            last_command_index = self.file_menu.index("end")
            if last_command_index == 0: # first label is deleted in search tab, create again.
                self.file_menu.insert_command(0, label=new_label, command=self.menu_edit_motif)
            else:
                self.file_menu.entryconfig(0, label=new_label)

    def get_menu_commands(self, menu):
        """
        Return commands under a menu in a list.
        :param menu:  menu is any level menu
        :return:
        """

        # Get all commands in the file menu
        commands = []
        for i in range(self.file_menu.index("end") + 1):
            label = menu.entrycget(i, "label")
            commands.append(label)
        return commands

    def on_tab_change(self, event):
        pass
        # selected_tab = event.widget.tab(event.widget.index("current"), "text")
        # print(f"Selected tab: {selected_tab}")
        #
        # # Perform any specific action based on the selected tab
        # if selected_tab == self.notebook_tabs[0]:
        #     self.change_menu_edit_button(self.menu_buttons[0][0])
        # elif selected_tab == self.notebook_tabs[1]:
        #     self.change_menu_edit_button(self.menu_buttons[1])
        # elif selected_tab == self.notebook_tabs[2]:
        #     self.change_menu_edit_button(self.menu_buttons[2])
        # elif selected_tab == self.notebook_tabs[3]:
        #     self.change_menu_edit_button(self.menu_buttons[3])
        # elif selected_tab == self.notebook_tabs[4]:
        #     self.change_menu_edit_button(None)

    def save_motif(self):
        pass

    def toggle_residue(self, *args):
        if self.motif_edit_on and not self.motif_edit_instance_mode:
            self.recreate_counts_table()

    def toggle_length_type(self, *args):
        if self.selected_length_type.get() == self.len_type_opts[0]:
            self.length_selector.config(state=tk.NORMAL)
        else:
            self.selected_length.set(0)
            self.length_selector.config(state=tk.DISABLED)

    def toggle_length(self):
        if self.motif_edit_on and not self.motif_edit_instance_mode:
            self.recreate_counts_table()

    def recreate_counts_table(self, from_scratch=False):
        residue = self.selected_residue.get()
        length = self.selected_length.get()

        # new motif
        alphabet = DNA_BASES if residue == self.residue_opts[0] else (
            RNA_BASES if residue == self.residue_opts[1] else AA_BASES)

        # remove old Counts table (if colns are not same), or empty (if colns are same)
        if not from_scratch and self.counts_table is not None:
            # Check column names are suitable
            cids = self.counts_table["columns"]
            column_names = [self.counts_table.heading(cid)["text"] for cid in cids]
            new_coln_names = ["Position"] + [letter for letter in alphabet]

            if column_names == new_coln_names:
                self.adjust_table_row_count(self.counts_table, length)

            elif len(column_names) == len(new_coln_names):
                self.rename_treeview_columns(self.counts_table, alphabet)  # rename columns
                self.adjust_table_row_count(self.counts_table, length)
            else:
                self.create_new_table_and_fill_with_zero(alphabet, length)
        else:
            self.create_new_table_and_fill_with_zero(alphabet, length)

        # pack
        self.counts_table.pack(expand=True, fill="both")

        # select first row
        self.select_first_row(self.counts_table)

    def adjust_table_row_count(self, tree, new_count):
        """
        Add new rows at the bottom or delete rows from bottom.
        :param tree:
        :param new_count:
        :return:
        """
        row_count = len(tree.get_children())
        if new_count < row_count:
            self.clear_treeview(tree, row_count-new_count)
        elif new_count > row_count:
            self.extend_treeview(tree, new_count-row_count)

    def create_new_table_and_fill_with_zero(self, alphabet, no_of_rows):
        # remove old treeview if exist
        if self.counts_table is not None:
            self.counts_table.destroy()

        # create new table and fill with zeros
        self.counts_table = self.create_treeview(self.counts_tab, ("Position", *alphabet))

        self.extend_treeview(self.counts_table, no_of_rows)

        # values_generator = ((i + 1, *pos_counts) for i, (*pos_counts,) in enumerate(zip(*counts.values())))
        # self.populate_treeview(self.counts_table, values_generator)

    def toggle_start_options(self):
        if self.search_scope_var.get() == "all":
            self.start_selected_radio.config(state=tk.NORMAL)
            self.start_beginning_radio.config(state=tk.NORMAL)
        else:
            self.start_selected_radio.config(state=tk.DISABLED)
            self.start_beginning_radio.config(state=tk.DISABLED)

    def toggle_threshold(self):
        default = 1
        if self.search_type_var.get() == "pssm":
            self.threshold_entry.config(state=tk.NORMAL)
            self.threshold_entry.delete(0, tk.END)
            self.threshold_entry.insert(0, str(default))
        else:
            self.threshold_entry.delete(0, tk.END)
            self.threshold_entry.config(state=tk.DISABLED)

    def disable_menu(self):
        toggle_menu_on_off(self.menu_bar, ["File"], False)

    def enable_menu(self):
        toggle_menu_on_off(self.menu_bar, ["File"], True)


class EditPseudoCounts(ModalWindow):
    def __init__(self, parent, curr_motif):
        self.parent = parent
        self.curr_motif = curr_motif
        super().__init__(self.parent, resizable=False)

        # Create a new window for editing
        self.title("PseudoCounts")

        # Add labels and entry widgets for each column
        self.entry_widgets = []
        row_offset = 0
        col_offset = 0
        for index, letter in enumerate(self.curr_motif.alphabet):
            if index == 10:  # for AA, the rest 10 letters are placed at right column
                row_offset = 10
                col_offset = 3

            label = tk.Label(self, text=letter)
            label.grid(row=index-row_offset, column=0 + col_offset, padx=10, pady=5)

            entry = Spinbox(self, from_=0, to=10000, width=5)
            entry.grid(row=index-row_offset, column=1 + col_offset, padx=20, pady=5)
            entry.insert(0, self.curr_motif.pseudocounts[letter])

            self.entry_widgets.append(entry)  # Store entry widgets in the new list

        # Add a "Save" button
        save_button = tk.Button(self, text="Save", command=self.save_changes)
        save_button.grid(row=len(self.curr_motif.alphabet), columnspan=2 + col_offset, pady=10)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    # Function to save the new values and update the Treeview
    def save_changes(self):
        new_values = [entry.get() for entry in self.entry_widgets]
        pseudo_counts = {letter: float(count) for letter, count in zip(self.curr_motif.alphabet, new_values)}
        print(pseudo_counts)
        self.curr_motif.pseudocounts = pseudo_counts
        self.destroy()


class EditBackground(ModalWindow):
    def __init__(self, parent, curr_motif):
        self.parent = parent
        self.curr_motif = curr_motif
        super().__init__(self.parent, resizable=False)

        # Create a new window for editing
        self.title("Background")

        # Add labels and entry widgets for each column
        self.entry_widgets = []
        row_offset = 0
        col_offset = 0
        for index, letter in enumerate(self.curr_motif.alphabet):
            if index == 10:  # for AA, the rest 10 letters are placed at right column
                row_offset = 10
                col_offset = 3

            label = tk.Label(self, text=letter)
            label.grid(row=index-row_offset, column=0 + col_offset, padx=10, pady=5)

            entry = tk.Entry(self, width=5)
            entry.grid(row=index-row_offset, column=1 + col_offset, padx=20, pady=5)
            entry.insert(0, self.curr_motif.background[letter])

            self.entry_widgets.append(entry)  # Store entry widgets in the new list

        # Add a "Save" button
        save_button = tk.Button(self, text="Save", command=self.save_changes)
        save_button.grid(row=len(self.curr_motif.alphabet), columnspan=2 + col_offset, pady=10)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    def disable_menu(self):
        pass

    def enable_menu(self):
        pass

    # Function to save the new values and update the Treeview
    def save_changes(self):
        # check if values are float
        new_values = [entry.get() for entry in self.entry_widgets]
        new_values_f = []
        for value in new_values:
            try:
                float_value = float(value)
                new_values_f.append(float_value)
            except ValueError:
                TKMessageBox("Error", f"Please enter a (decimal) number.", messagebox.showerror,
                             self).show()
                return

        # check sign
        if any([value < 0 for value in new_values_f]):
            TKMessageBox("Error", f"Please enter a non - negative number.", messagebox.showerror,
                         self).show()
            return

        # check sum
        if sum(new_values_f) != 1:
            TKMessageBox("Error", f"Given numbers must add up to 1.", messagebox.showerror,
                         self).show()
            return

        background = {letter: bg for letter, bg in zip(self.curr_motif.alphabet, new_values_f)}
        print(background)
        self.curr_motif.background = background
        self.destroy()


class EditSelectedRow(ModalWindow):
    def __init__(self, parent, tree, disabled_colns):
        self.parent = parent
        self.tree = tree
        self.disabled_colns = disabled_colns
        super().__init__(self.parent, resizable=False)

        self.selected_item = self.tree.focus()  # Get the ID of the selected item
        if not self.selected_item:
            self.destroy()
        selected_values = self.tree.item(self.selected_item, 'values')
        self.position = selected_values[0]

        # Create a new window for editing
        self.title("Edit Row")

        # Add labels and entry widgets for each column
        self.entry_widgets = []
        row_offset = 0
        col_offset = 0
        column_names = [self.tree.heading(cid)["text"] for cid in self.tree["columns"]]
        for col_index, col_name in enumerate(column_names[1:]):
            if col_index == 10:  # for AA, the rest 10 letters are placed at right column
                row_offset = 10
                col_offset = 3

            label = tk.Label(self, text=col_name)
            label.grid(row=col_index-row_offset, column=0 + col_offset, padx=5, pady=5)

            entry = tk.Entry(self)
            entry.grid(row=col_index-row_offset, column=1 + col_offset, padx=5, pady=5)
            entry.insert(0, selected_values[col_index+1])
            if col_name in self.disabled_colns:  # unnecessary
                entry.config(state="readonly")
            self.entry_widgets.append(entry)  # Store entry widgets in the new list

        # Add a "Save" button
        save_button = tk.Button(self, text="Save", command=self.save_changes)
        save_button.grid(row=len(tree["columns"]), columnspan=2 + col_offset, pady=10)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    def disable_menu(self):
        pass

    def enable_menu(self):
        pass

    # Function to save the new values and update the Treeview
    def save_changes(self):

        # check if values are float
        new_values = [entry.get() for entry in self.entry_widgets]
        new_values_f = []
        for value in new_values:
            try:
                float_value = float(value)
                new_values_f.append(float_value)
            except ValueError:
                TKMessageBox("Error", f"Please enter a (decimal) number.", messagebox.showerror,
                             self).show()
                return

        # check sign
        if any([value < 0 for value in new_values_f]):
            TKMessageBox("Error", f"Please enter a non - negative number.", messagebox.showerror,
                         self).show()
            return


        new_values = [self.position]+new_values_f
        self.tree.item(self.selected_item, values=new_values)
        #print_all_table_items(self.tree)
        self.destroy()


class SelectEditMode(ModalWindow):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(self.parent, resizable=False)

        self.opts = ["Instances", "Counts"]
        self.is_instance_mode = True
        self.title("Select Edit Mode")
        self.is_closed_succesfully = False

        self.option_var = tk.StringVar(value=self.opts[0])

        # Create radio buttons for the options
        self.instances_radio = ttk.Radiobutton(self, text=self.opts[0], variable=self.option_var, value=self.opts[0])
        self.counts_radio = ttk.Radiobutton(self, text=self.opts[1], variable=self.option_var, value=self.opts[1])

        self.instances_radio.pack(padx = 50, pady=5)
        self.counts_radio.pack(padx = 50, pady=5)

        # Create a button to confirm the selection
        self.confirm_button = ttk.Button(self, text="Confirm", command=self.confirm_selection)
        self.confirm_button.pack(pady=10)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    def confirm_selection(self):
        self.is_instance_mode = self.option_var.get() == self.opts[0]
        self.is_closed_succesfully = True
        self.destroy()


class SaveMotifDialog(ModalWindow):
    def __init__(self, parent, motifs):
        self.parent = parent
        self.motifs = motifs
        super().__init__(self.parent, resizable=False)
        self.title("Select Items")

        self.selected_motifs = None
        self.selected_format = None
        self.selected_path = None
        self.is_closed_succesfully = False

        self.supported_formats = [
            "clusterbuster",
            "pfm",
            "jaspar",
            "transfac"
        ]

        top_frame = tk.Frame(self)
        top_frame.pack(fill="both", expand=True, pady=10)

        # Create a Canvas widget to enable scrolling
        canvas = tk.Canvas(top_frame)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Create a vertical scrollbar linked to the Canvas
        scrollbar = tk.Scrollbar(top_frame, orient=tk.VERTICAL, command=canvas.yview)
        scrollbar.place(relx=1, rely=0, relheight=1, anchor='ne')

        # Link the canvas and the scrollbar
        canvas.config(yscrollcommand=scrollbar.set)

        # Create a frame to contain the checkbuttons
        checkbutton_frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=checkbutton_frame, anchor="nw")

        checkbutton_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        # Populate the frame with many Checkbuttons
        self.checkbuttons = {}
        for motif_id, motif in self.motifs.items():
            var = tk.BooleanVar(value=True)
            residue_type = "DNA" if motif.alphabet == DNA_BASES else ("RNA" if motif.alphabet == RNA_BASES else "Aminoacid")
            checkbutton = tk.Checkbutton(checkbutton_frame, text=motif_id + f" ({residue_type})", variable=var)
            checkbutton.pack(anchor='w')
            self.checkbuttons[motif_id] = var

        # Bottom frame for path and format
        bottom_frame = tk.Frame(self)
        bottom_frame.pack(pady=10)

        # bottom_frame.columnconfigure(0, weight=1, uniform="a")
        # bottom_frame.columnconfigure(1, weight=2, uniform="a")
        # bottom_frame.columnconfigure(2, weight=1, uniform="a")
        # bottom_frame.rowconfigure(0, weight=1)
        # bottom_frame.rowconfigure(1, weight=1)
        # bottom_frame.rowconfigure(2, weight=1)

        tk.Label(bottom_frame, text="Path:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.output_file_entry = tk.Entry(bottom_frame, width=20)
        self.output_file_entry.grid(row=0, column=1, padx=5, pady=5, sticky='we')
        browse = tk.Button(bottom_frame, text="Browse", command=self.browse_output_file)
        browse.grid(row=0, column=2, padx=5, pady=5)

        tk.Label(bottom_frame, text="Format:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.format_var = tk.StringVar(self)
        self.format_var.set(self.supported_formats[0])
        self.format_menu = tk.OptionMenu(bottom_frame, self.format_var, *self.supported_formats)
        self.format_menu.grid(row=1, column=1, padx=5, pady=5, sticky='w')

        save_button = tk.Button(bottom_frame, text="Save", command=self.save_selected_motifs)
        save_button.grid(row=2, column=1, padx=5, pady=5)

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

    def save_selected_motifs(self):
        """
        Notes:
        - "clusterbuster" does not support aminoacid and RNA.
        - "pfm" does not support aminoacid and RNA. It does not support multiple motifs.
        - "jaspar" does not support aminoacid and RNA.
        - "transfac" supports everything.

        :return:
        """
        selected_motif_ids = [motif_id for motif_id, var in self.checkbuttons.items() if var.get()]
        if len(selected_motif_ids) == 0:
            TKMessageBox("Error", f"Please select at least one motif.", messagebox.showerror, self).show()
            return

        self.selected_motifs = [self.motifs[motif_id] for motif_id in selected_motif_ids]
        self.selected_format = self.format_var.get()
        self.selected_path = self.output_file_entry.get()

        directory = os.path.dirname(self.selected_path)
        if not os.path.exists(directory):
            TKMessageBox("Error", f"Please select a valid path.", messagebox.showerror, self).show()
            return

        if self.selected_format == self.supported_formats[1] and len(self.selected_motifs) > 1:
            TKMessageBox(
                "Error",
                f"{self.selected_format.upper()} format does not support multiple motifs. Please select a single motif or another format.",
                messagebox.showerror,
                self
            ).show()
            return

        if self.selected_format in self.supported_formats[:3] and any(motif.alphabet != DNA_BASES for motif in self.selected_motifs) :
            TKMessageBox(
                "Error",
                f"{self.selected_format.upper()} format only supports DNA motifs.",
                messagebox.showerror,
                self
            ).show()
            return

        self.is_closed_succesfully = True
        self.destroy()


class SearchManager:
    def __init__(self, seq_viewer: SequenceViewer, motif_finder: MotifFinder ,curr_search_motif, pssm_search: bool, pssm_threshold: float, search_data:list, start_line_index: int, stop_line_index: int):
        self.seq_viewer = seq_viewer
        self.motif_finder = motif_finder
        self.curr_search_motif = curr_search_motif
        self.pssm_search = pssm_search
        self.pssm_threshold = pssm_threshold
        self.search_data = search_data
        self.start_line_index = start_line_index
        self.stop_line_index = stop_line_index

        # if true, search completed, user can browse among found matches.
        # if false, generator still continues to find new.
        self.search_completed = False

        # this changes only with "start" and "stop" buttons
        self.search_is_on = True

        # query thread shared variable
        self.query_thread = None
        self.query_result = None  # this is reset to None after finding each match

        if isinstance(self.curr_search_motif, Motif):
            self.fixed_length = True
            self.length = len(self.curr_search_motif)
        else:
            self.fixed_length = False
            self.length = None

        # init
        self.search_curr_seq_index = 0  # the sequence which is being searched currently.
        self.search_matches = []  # (line_number, (pos, instance/pssm score))
        self.curr_match_index = -1  # currently shown match

        self.curr_search_gen = self.create_generator(self.search_data[0])

    def create_generator(self, line):  # line: Seq or str
        line = line.upper()
        if self.pssm_search:  # PSSM
            curr_search_gen = self.curr_search_motif.pssm.search(line, threshold=self.pssm_threshold, both=False)

        else:  # EXACT
            if self.fixed_length:  # fixed length
                curr_search_gen = self.curr_search_motif.instances.search(line)

            else:  # variable length
                s = Seq(line) if not isinstance(line, Seq) else line
                curr_search_gen = s.search(self.curr_search_motif.instances)
        return curr_search_gen

    def prev(self):
        self.curr_match_index -= 1
        match = self.search_matches[self.curr_match_index]
        print(self.curr_match_index, match)
        return match

    def next(self):
        """
        Searches until no match left in all search region. If found, return first match.
        If already found match is demanded (due to some "prev" clicks), return from list. If not, run a thread to find new.
        :return: None if no left.
        """
        # query already found match (due to some "prev" clicks)
        if len(self.search_matches) >= self.curr_match_index + 2:
            self.curr_match_index += 1
            match = self.search_matches[self.curr_match_index]
            print(self.curr_match_index, match)
            return match

        # Start new thread
        self.query_result = None
        self.query_thread = threading.Thread(target=self.query_new_match)
        self.query_thread.start()
        return None  # starts new thread and returns None

    def query_new_match(self):
        """
        Runs a query for next item in a new thread. After completion, thread order "process_query_result" or
        "stop_searching_gui" functions to be run by main thread based on condition.
        based on
        :return:
        """
        # query new match
        while True:
            # check if user clicked "stop" button
            with self.motif_finder.stop_search_lock:
                if not self.search_is_on:
                    break  # user clicked "stop"

            try:
                new_found = next(self.curr_search_gen)  # return next match (triggers biopython generator)
                new_match = self.create_new_match(new_found)
                self.search_matches.append(new_match)
                self.curr_match_index += 1
                print(self.curr_match_index, new_match)
                self.query_result = new_match
                self.motif_finder.after(0, self.motif_finder.process_query_result)  # call on main thread
                return  # We found new match, stop thread.

            except StopIteration:
                self.search_curr_seq_index += 1  # proceed to next seq
                if self.search_curr_seq_index <= len(self.search_data) - 1:
                    print("New line:", self.search_curr_seq_index)
                    line = self.search_data[self.search_curr_seq_index]
                    self.curr_search_gen = self.create_generator(line)
                else:
                    break  # no line left

        # User stopped search using "STOP" button while searching
        with self.motif_finder.stop_search_lock:
            if not self.search_is_on:
                print("Stopped by user while searching.")
                self.query_result = None
                self.motif_finder.after(0, self.motif_finder.stop_searching_gui)  # call on main thread
                return
            # Search completed normally but NO MATCH left.
            else:
                self.search_completed = True
                print("Search completed. No next item exist.")
                self.query_result = None
                self.motif_finder.after(0, self.motif_finder.process_query_result)  # call on main thread
                return

    def create_new_match(self, new_found):
        if self.pssm_search:
            pos, pssm_score = new_found[0], new_found[1]
            line_index = self.search_curr_seq_index + self.start_line_index  # note that we add offset
            new_match = self.Match(line_index, pos, length=self.length, pssm_score=pssm_score)
        else:
            pos, instance = new_found[0], new_found[1]
            line_index = self.search_curr_seq_index + self.start_line_index  # note that we add offset
            new_match = self.Match(line_index, pos, length=len(instance), instance=instance)
        return new_match

    def get_match_count(self):
        return len(self.search_matches) if self.search_completed else None

    class Match:
        def __init__(self, line_index, pos, length, instance=None, pssm_score=None):
            if all(var is None for var in [instance, pssm_score]):
                raise Exception("Give instance or pssm_score.")

            if all(var is not None for var in [instance, pssm_score]):
                raise Exception("Give instance or pssm_score.")

            self.line_index = line_index
            self.pos = pos
            self.length = length
            self.instance = instance
            self.pssm_score = pssm_score

        def __str__(self):
            return (f"line_index: {self.line_index}, "
                    f"pos: {self.pos}, "
                    f"length: {self.length}, "
                    f"instance: {self.instance}, "
                    f"pssm_score: {self.pssm_score}")