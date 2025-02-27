import tkinter as tk
from tkinter import ttk, messagebox
from Bio.SeqFeature import SeqFeature, Reference, Location, CompoundLocation, SimpleLocation, UnknownPosition, \
    UncertainPosition, OneOfPosition, WithinPosition, AfterPosition, BetweenPosition, BeforePosition, ExactPosition
from utils.model_window import ModalWindow
from utils.utils import toggle_menu_on_off
from utils.utils import print_all_table_items
from utils.windows import TKMessageBox
from typing import Callable
import copy

special_field_types = {
    "molecule_type": str,
    "data_file_division": str,
    "date": str,
    "accessions": list,
    "sequence_version": int,
    "gi": str,
    "keywords": list,
    "source": str,
    "organism": str,
    "taxonomy": list,
    "references": list,
    "comment": str
}

# Mapping for position hints
POSITION_HINTS = {
    "ExactPosition": "pos",
    "BeforePosition": "pos",
    "AfterPosition": "pos",
    "WithinPosition": "start_pos,end_pos",
    "BetweenPosition": "start_pos,end_pos",
    "OneOfPosition": "pos1,pos2,pos3,...",
    "UncertainPosition": "pos",
    "UnknownPosition": ""
}

POSITION_VALUES = {
    "ExactPosition": "Exact Position",
    "BeforePosition": "Before Position",
    "AfterPosition": "After Position",
    "WithinPosition": "Within Position",
    "BetweenPosition": "Between Position",
    "OneOfPosition": "One Of Position",
    "UncertainPosition": "Uncertain Position",
    "UnknownPosition": "Unknown Position"
}


class AnnotationViewer(ModalWindow):
    def __init__(self, seq_id, seq_record, line_text, parent, show_letter_annot: bool, edit_callback:Callable, save_callback:Callable, line_no):
        self.parent = parent  # parent window
        super().__init__(self.parent, min_size=(400, 300))
        self.seq_id = seq_id
        self.seq_record = seq_record

        # seq annot data
        self.line_text = line_text  # seq in tw (different that seqrecord if edited)
        self.id = self.seq_record.id
        self.name = self.seq_record.name
        self.description = self.seq_record.description
        self.annotations = copy.deepcopy(self.seq_record.annotations)  # apply deepcopy as we modify it
        self.features = copy.deepcopy(self.seq_record.features)  # apply deepcopy as we modify it
        self.compound_loc_to_list()
        self.letter_annots = self.seq_record.letter_annotations if show_letter_annot else dict()  # Empty if seq modified (ONLY LETTER ANNOT EMPTIED)
        self.letter_annots_types = {key: type(value) for key, value in self.letter_annots.items()}
        self.dbxrefs = self.seq_record.dbxrefs

        self.edit_callback = edit_callback
        self.save_callback = save_callback
        self.line_no = line_no

        # Configure window
        self.title("Annotation Viewer")

        # Create menu
        self.menu_bar = tk.Menu(self)
        self.menu_names = ["File", "Edit"]
        self.file_menu = None
        self.edit_menu = None
        self.create_menu()

        # Create tables in Notebook
        self.notebook = ttk.Notebook(self)
        self.tree_info, self.tree_nov_annot, self.tree_dbx_refs = None, None, None # normal trees
        self.tree_nov_feature, self.tree_letter_annot = None, None  # nested trees
        self.create_tables()
        self.notebook.pack(fill='both', expand=True)

        # Bind the <<NotebookTabChanged>> event
        self.notebook.bind('<<NotebookTabChanged>>', self.on_tab_change)
        self.selected_tab = None

        self.tab_names = ["Information", "Annotation", "Features", "Letter Annotation", "Database cross-references"]

        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def on_close(self):
        quit = TKMessageBox("Warning", "Are you sure to quit now ?  Unsaved annots will be lost.", messagebox.askyesno,
                            self).show()
        if quit:
            self.destroy()
        else:
            return  # kapatmaktan vazgeçti

    def compound_loc_to_list(self):  # convert CompoundLocation to list of SimpleLocation
        for feature in self.features:
            if isinstance(feature.location, CompoundLocation):
                feature.location = [simple_loc for simple_loc in feature.location.parts]

        if 'references' in self.annotations:
            for ref in self.annotations['references']:
                for i, loc in enumerate(ref.location):
                    if isinstance(loc, CompoundLocation):
                        ref.location[i] = [simple_loc for simple_loc in loc.parts]

    def list_to_compound_loc(self):  # convert list of SimpleLocation to CompoundLocation
        for feature in self.features:
            if isinstance(feature.location, list):
                feature.location = CompoundLocation(parts=feature.location)

        if 'references' in self.annotations:
            for ref in self.annotations['references']:
                for i, loc in enumerate(ref.location):
                    if isinstance(loc, list):
                        ref.location[i] = CompoundLocation(parts=loc)

    def create_menu(self):
        # File Menu
        self.file_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.file_menu.add_command(label="Save and Close", command=self.save_and_close)
        self.file_menu.add_command(label="Discard and Close", command=self.discard_and_close)

        # Edit Menu
        self.edit_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.edit_menu.add_command(label="Edit selected row", command=self.edit_selected_row)
        self.edit_menu.add_command(label="Add new row", command=self.add_new_row)
        self.edit_menu.add_command(label="Add sub-row", command=self.add_sub_row)

        self.genbank_menu = tk.Menu(self.edit_menu, tearoff=0)
        for field in special_field_types:
            self.genbank_menu.add_command(label=field, command=lambda x=field: self.add_special_row(x))
        self.edit_menu.add_cascade(label="Add Genbank-EBML row", menu=self.genbank_menu)

        self.edit_menu.add_command(label="Delete selected row", command=self.delete_selected_row)
        # self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Add new column", command=self.add_new_column)
        self.edit_menu.add_command(label="Delete column", command=self.delete_column)

        # Add menus to the menu bar
        self.menu_bar.add_cascade(label=self.menu_names[0], menu=self.file_menu)
        self.menu_bar.add_cascade(label=self.menu_names[1], menu=self.edit_menu)

        # Configure the root window with the menu bar
        self.config(menu=self.menu_bar)

    def disable_menu(self):
        toggle_menu_on_off(self.menu_bar, self.menu_names, False)

    def enable_menu(self):
        toggle_menu_on_off(self.menu_bar, self.menu_names, True)

    def create_tables(self):
        # tab 1 - Info table
        tab_info = ttk.Frame(self.notebook)
        self.tree_info = self.create_treeview(tab_info, ("Field", "Value"))
        self.populate_treeview(self.tree_info, [("Id", self.id),
                                                ("Name", self.name),
                                                ("Description", self.description)])
        # self.tree_info.pack(expand=True, fill="both")
        self.tree_info.pack(side="left", fill="both", expand=True)
        self.notebook.add(tab_info, text='Information')

        # tab 2 - Annot table
        tab_annot = ttk.Frame(self.notebook)
        self.tree_nov_annot = NestedObjectViewer(self, tab_annot, self.annotations)
        self.tree_nov_annot.tree.pack(side="left", fill="both", expand=True)
        self.notebook.add(tab_annot, text='Annotation')

        # tab 3 - Feature table
        tab_feature = ttk.Frame(self.notebook)
        self.tree_nov_feature = NestedObjectViewer(self, tab_feature, self.features)
        self.tree_nov_feature.tree.pack(side="left", fill="both", expand=True)
        self.notebook.add(tab_feature, text='Features')

        # tab 4 - Letter annot table
        tab_letter = ttk.Frame(self.notebook)
        self.tree_letter_annot = self.create_treeview(tab_letter,
                                                      ("Index", "Letter", *self.letter_annots.keys()))
        values_generator = ((i + 1, seq_char, *annotations) for i, (seq_char, *annotations) in
                            enumerate(zip(self.line_text, *self.letter_annots.values())))
        self.populate_treeview(self.tree_letter_annot, values_generator)
        # self.tree_letter_annot.pack(expand=True, fill="both")
        self.tree_letter_annot.pack(side="left", fill="both", expand=True)
        self.notebook.add(tab_letter, text='Letter Annotation')

        # tab 5 - DBX refs table
        tab_dbx = ttk.Frame(self.notebook)
        self.tree_dbx_refs = self.create_treeview(tab_dbx, ("Value"))
        self.populate_treeview(self.tree_dbx_refs, self.dbxrefs)
        # self.tree_dbx_refs.pack(expand=True, fill="both")
        self.tree_dbx_refs.pack(side="left", fill="both", expand=True)
        self.notebook.add(tab_dbx, text='Database cross-references')

        # # Bind the key press event
        # self.tree.bind('<KeyPress>', lambda event: self.on_key_pressed(event))
        #
        # # Pack the Treeview
        # self.tree.pack(expand=True, fill="both")

    def on_key_pressed(self, event):
        # Check if the pressed key is "e" (case-insensitive) and if a row is selected
        if event.char.lower() == "e":
            selected_item = self.tree_info.selection()
            if selected_item:
                print(f"Key 'E' pressed after selecting row: {selected_item[0]}")

    def create_treeview(self, tab, columns):
        tree = ttk.Treeview(tab, show='headings', selectmode='browse', columns=columns)
        for column in list(tree['columns']):
            tree.heading(column, text=column)
            tree.column(column, minwidth=100)

        # Create vertical scrollbar
        vscrollbar = ttk.Scrollbar(tab, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=vscrollbar.set)
        # vscrollbar.place(relx=1, rely=0, relheight=1, anchor='ne')
        vscrollbar.pack(side="right", fill="y")

        # Create horizontal scrollbar
        hscrollbar = ttk.Scrollbar(tab, orient="horizontal", command=tree.xview)
        tree.configure(xscrollcommand=hscrollbar.set)
        # hscrollbar.place(relx=0, rely=1, relwidth=1, anchor='sw')
        hscrollbar.pack(side="bottom", fill="x")

        # Bind the Configure event to adjust column widths
        tree.bind("<Configure>", lambda event: self.tree_on_configure(event, tree))

        return tree

    def populate_treeview(self, tree, data_iterator):
        for values in data_iterator:
            tree.insert("", index="end", values=values)

    def tree_on_configure(self, event, tree):
        # Adjust column widths when the window size changes
        for column in tree['columns']:
            tree.column(column, width=200, minwidth=100)

    def check_selected(self, tree):
        if len(tree.selection()) == 0:
            TKMessageBox("Error", "No selected row !", messagebox.showerror, self).show()
            return False
        return True

    def edit_selected_row(self):
        if self.selected_tab == self.tab_names[0]:
            print("Edit selected row (Information Menu)")

            if self.check_selected(self.tree_info):
                self.selected_item = self.tree_info.selection()[0]  # Get the ID of the selected item
                selected_values = self.tree_info.item(self.selected_item, 'values')
                if selected_values[0] == "Description":
                    EditSelectedRowLarge(self, self.tree_info)
                else:
                    EditSelectedRow(self, self.tree_info, ["Field"])

        elif self.selected_tab == self.tab_names[1]:
            print("Edit selected row (Annotations Menu)")

            if self.check_selected(self.tree_nov_annot.tree):
                self.tree_nov_annot.edit_selected_row()

        elif self.selected_tab == self.tab_names[2]:
            print("Edit selected row (Features Menu)")

            if self.check_selected(self.tree_nov_feature.tree):
                self.tree_nov_feature.edit_selected_row()

        elif self.selected_tab == self.tab_names[3]:
            print("Edit selected row (Letter Annotation Menu)")

            if self.check_selected(self.tree_letter_annot):
                EditSelectedRowLetterAnnot(self, self.tree_letter_annot, ["Index"])

        elif self.selected_tab == self.tab_names[4]:
            print("Edit selected row (Database Cross-references Menu)")

            if self.check_selected(self.tree_dbx_refs):
                EditSelectedRow(self, self.tree_dbx_refs, [])

    def add_new_row(self):
        if self.selected_tab == self.tab_names[1]:
            print("Add new row (Annotations Menu)")
            self.tree_nov_annot.add_new_row(True)

        elif self.selected_tab == self.tab_names[2]:
            print("Add new row (Features Menu)")
            self.tree_nov_feature.add_new_row(False)

        elif self.selected_tab == self.tab_names[3]:
            print("Add new row (Letter Annotations Menu)")
            AddNewRowLetterAnnot(self, self.tree_letter_annot)

        elif self.selected_tab == self.tab_names[4]:
            print("Add new row (Database Cross-references Menu)")
            AddNewRow(self, self.tree_dbx_refs)

    def add_sub_row(self):
        if self.selected_tab == self.tab_names[1]:
            print("Add sub row (Annotations Menu)")
            self.tree_nov_annot.add_sub_row(True)

        elif self.selected_tab == self.tab_names[2]:
            print("Add sub row (Features Menu)")
            self.tree_nov_feature.add_sub_row(False)

    def add_special_row(self, field):
        self.tree_nov_annot.add_special_row(field)

    def delete_selected_row(self):
        if self.selected_tab == self.tab_names[1]:
            print("Delete selected row (Annotations Menu)")

            if self.check_selected(self.tree_nov_annot.tree):
                self.tree_nov_annot.delete_selected_row()

        elif self.selected_tab == self.tab_names[2]:
            print("Delete selected row (Features Menu)")

            if self.check_selected(self.tree_nov_feature.tree):
                self.tree_nov_feature.delete_selected_row()

        elif self.selected_tab == self.tab_names[3]:
            print("Delete selected row (Letter Annotations Menu)")

            if self.check_selected(self.tree_letter_annot):
                delete_selected_row(self.tree_letter_annot)
                update_indices(self.tree_letter_annot)

        elif self.selected_tab == self.tab_names[4]:
            print("Delete selected row (Database Cross-references Menu)")

            if self.check_selected(self.tree_dbx_refs):
                delete_selected_row(self.tree_dbx_refs)

    def add_new_column(self):
        if self.selected_tab == self.tab_names[3]:
            print("Add new column (Letter Annotation Menu)")
            AddNewColumn(self, self.tree_letter_annot)

    def delete_column(self):
        if self.selected_tab == self.tab_names[3]:
            print("Delete column (Letter Annotation Menu)")
            tree = self.tree_letter_annot
            current_columns = [tree.heading(key)["text"] for key in list(tree['columns'][2:])]  # get heading info
            if len(current_columns) > 0:
                DeleteColumn(self, self.tree_letter_annot, current_columns)
            else:
                TKMessageBox("Error", "No column to delete", messagebox.showerror, self).show()

    def on_tab_change(self, event):
        self.selected_tab = event.widget.tab(event.widget.index("current"), "text")

        # enable-disable edit menu buttons
        if self.selected_tab == self.tab_names[0]:  # info
            self.edit_menu.entryconfig(0, state=tk.NORMAL)  # Edit selected row
            self.edit_menu.entryconfig(1, state=tk.DISABLED)  # add new row
            self.edit_menu.entryconfig(2, state=tk.DISABLED)  # add sub row
            self.edit_menu.entryconfig(3, state=tk.DISABLED)  # add special row
            self.edit_menu.entryconfig(4, state=tk.DISABLED)  # delete selected row
            self.edit_menu.entryconfig(5, state=tk.DISABLED)  # add new column
            self.edit_menu.entryconfig(6, state=tk.DISABLED)  # delete column
        elif self.selected_tab == self.tab_names[1]:  # annot
            self.edit_menu.entryconfig(0, state=tk.NORMAL)  # Edit selected row
            self.edit_menu.entryconfig(1, state=tk.NORMAL)  # add new row
            self.edit_menu.entryconfig(2, state=tk.NORMAL)  # add sub row
            self.edit_menu.entryconfig(3, state=tk.NORMAL)  # add special row
            self.edit_menu.entryconfig(4, state=tk.NORMAL)  # delete selected row
            self.edit_menu.entryconfig(5, state=tk.DISABLED)  # add new column
            self.edit_menu.entryconfig(6, state=tk.DISABLED)  # delete column
        elif self.selected_tab == self.tab_names[2]:  # features
            self.edit_menu.entryconfig(0, state=tk.NORMAL)  # Edit selected row
            self.edit_menu.entryconfig(1, state=tk.NORMAL)  # add new row
            self.edit_menu.entryconfig(2, state=tk.NORMAL)  # add sub row
            self.edit_menu.entryconfig(3, state=tk.DISABLED)  # add special row
            self.edit_menu.entryconfig(4, state=tk.NORMAL)  # delete selected row
            self.edit_menu.entryconfig(5, state=tk.DISABLED)  # add new column
            self.edit_menu.entryconfig(6, state=tk.DISABLED)  # delete column
        elif self.selected_tab == self.tab_names[3]:  # letter annot
            self.edit_menu.entryconfig(0, state=tk.NORMAL)  # Edit selected row
            self.edit_menu.entryconfig(1, state=tk.NORMAL)  # add new row
            self.edit_menu.entryconfig(2, state=tk.DISABLED)  # add sub row
            self.edit_menu.entryconfig(3, state=tk.DISABLED)  # add special row
            self.edit_menu.entryconfig(4, state=tk.NORMAL)  # delete selected row
            self.edit_menu.entryconfig(5, state=tk.NORMAL)  # add new column
            self.edit_menu.entryconfig(6, state=tk.NORMAL)  # delete column
        elif self.selected_tab == self.tab_names[4]:  # database crx
            self.edit_menu.entryconfig(0, state=tk.NORMAL)  # Edit selected row
            self.edit_menu.entryconfig(1, state=tk.NORMAL)  # add new row
            self.edit_menu.entryconfig(2, state=tk.DISABLED)  # add sub row
            self.edit_menu.entryconfig(3, state=tk.DISABLED)  # add special row
            self.edit_menu.entryconfig(4, state=tk.NORMAL)  # delete selected row
            self.edit_menu.entryconfig(5, state=tk.DISABLED)  # add new column
            self.edit_menu.entryconfig(6, state=tk.DISABLED)  # delete column

    def save_and_close(self):
        """
        Read the tables and fills all annots of seq record with these values.
        Note: Always assign new attributes rather than modifying attributes internally.
        :return:
        """

        # Update Information tab
        info_values = [(self.tree_info.item(item_id, 'values')[0], self.tree_info.item(item_id, 'values')[1]) for item_id in self.tree_info.get_children()]
        for key, value in info_values:
            if key == "Id":
                self.seq_record.id = value
            elif key == "Name":
                self.seq_record.name = value
            elif key == "Description":
                self.seq_record.description = value

        # Update Annotation and Feature tab
        self.list_to_compound_loc()
        self.seq_record.annotations = self.annotations
        self.seq_record.features = self.features

        # Update seq and Letter Annotation attribute
        coln_offset = 2
        letters = [self.tree_letter_annot.item(item_id, 'values')[1] for item_id in self.tree_letter_annot.get_children()]
        letter_annot_values = [self.tree_letter_annot.item(item_id, 'values')[coln_offset:] for item_id in
                               self.tree_letter_annot.get_children()]
        letter_annot_dict = {key: values for key, values in
                             zip(self.tree_letter_annot["columns"][coln_offset:], zip(*letter_annot_values))}
        letter_annot_dict = self.convert_dict_values_to_float(letter_annot_dict)  # convert to float if possible
        self.seq_record.seq = "".join(letters)

        for key, value in letter_annot_dict.items():  # convert list to str if it is str initially (e.g. secondary structure)
            if key in self.letter_annots_types and self.letter_annots_types[key] == str:
                new_value = "".join(value)
                letter_annot_dict[key] = new_value

        self.seq_record.letter_annotations = letter_annot_dict

        # Update Database cross-references tab
        dbx_values = [self.tree_dbx_refs.item(item_id, 'values')[0] for item_id in self.tree_dbx_refs.get_children()]
        self.seq_record.dbxrefs = dbx_values

        # Edit text widget (seq id and seq content)
        self.edit_callback(self.line_no, self.seq_id, self.seq_record.id, self.seq_record.seq)

        # Call save callback
        self.save_callback()
        self.parent.unsaved_file = True
        self.destroy()

    def convert_dict_values_to_float(self, input_dict):
        try:
            converted_dict = {}
            for key, str_list in input_dict.items():
                converted_dict[key] = [float(item) for item in str_list]
        except ValueError:
            return input_dict  # original
        else:
            return converted_dict  # int

    def discard_and_close(self):
        self.destroy()


class EditSelectedRow(ModalWindow):
    def __init__(self, parent, tree, disabled_colns):
        self.parent = parent
        self.tree = tree
        self.disabled_colns = disabled_colns
        super().__init__(self.parent, resizable=(True, False), min_size=(200, 100))

        self.selected_item = self.tree.selection()[0]  # Get the ID of the selected item
        selected_values = self.tree.item(self.selected_item, 'values')

        # Create a new window for editing
        self.title("Edit Row")

        # Configure grid layout to allow resizing
        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)

        # Add labels and entry widgets for each column
        self.entry_widgets = []
        for col_index, col_name in enumerate(self.tree["columns"]):
            label = tk.Label(self, text=col_name)
            label.grid(row=col_index, column=0, padx=5, pady=5, sticky="w")

            # entry_var = tk.StringVar()  çalışmadı, get boş dünüyor.
            entry = tk.Entry(self)
            entry.grid(row=col_index, column=1, padx=5, pady=5, sticky="we")
            entry.insert(0, selected_values[col_index])
            if col_name in self.disabled_colns:
                entry.config(state="readonly")
            self.entry_widgets.append(entry)  # Store entry widgets in the new list

        # Add a "Save" button
        save_button = tk.Button(self, text="Save", command=self.save_changes)
        save_button.grid(row=len(tree["columns"]), columnspan=2, pady=10)

    def disable_menu(self):
        pass  # no menu exist

    def enable_menu(self):
        pass  # no menu exist

    # Function to save the new values and update the Treeview
    def save_changes(self):
        new_values = [entry.get() for entry in self.entry_widgets]
        self.tree.item(self.selected_item, values=new_values)
        print_all_table_items(self.tree)
        self.destroy()


class EditSelectedRowLarge(ModalWindow):
    def __init__(self, parent, tree):
        self.parent = parent
        self.tree = tree
        super().__init__(self.parent, min_size=(400, 200))

        self.selected_item = self.tree.selection()[0]  # Get the ID of the selected item
        selected_values = self.tree.item(self.selected_item, 'values')

        # Create a new window for editing
        self.title("Edit Row")
        self.geometry("400x200")

        # Configure grid layout to allow resizing
        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)

        # Field label + Entry
        ttk.Label(self, text="Field:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.field_entry = ttk.Entry(self)
        self.field_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=5)
        self.field_entry.insert(0, selected_values[0])
        self.field_entry.config(state="readonly")

        # Value label + Text
        ttk.Label(self, text="Value:").grid(row=1, column=0, sticky="nw", padx=5, pady=5)
        self.value_text = tk.Text(self, height=5, wrap="word")
        self.value_text.grid(row=1, column=1, sticky="nsew", padx=5, pady=5)
        self.value_text.insert("1.0", selected_values[1])

        # Add a "Save" button
        save_button = ttk.Button(self, text="Save", command=self.save_changes)
        save_button.grid(row=2,column=0, columnspan=2, pady=10)

        # Grid weight for resizing
        self.columnconfigure(1, weight=1)  # Allow horizontal resizing
        self.rowconfigure(1, weight=1)  # Allow vertical resizing for Text

    def disable_menu(self):
        pass  # no menu exist

    def enable_menu(self):
        pass  # no menu exist

    # Function to save the new values and update the Treeview
    def save_changes(self):
        new_values = [self.field_entry.get(), self.value_text.get("1.0", "end").strip()]
        self.tree.item(self.selected_item, values=new_values)
        print_all_table_items(self.tree)
        self.destroy()


class EditSelectedRowLetterAnnot(ModalWindow):
    def __init__(self, parent, tree, disabled_colns):
        self.parent = parent
        self.tree = tree
        self.disabled_colns = disabled_colns
        super().__init__(self.parent, resizable=False)

        self.selected_item = self.tree.selection()[0]  # Get the ID of the selected item
        selected_values = self.tree.item(self.selected_item, 'values')

        # Create a new window for editing
        self.title("Edit Row")

        # Add labels and entry widgets for each column
        self.entry_widgets = []
        for col_index, col_name in enumerate(self.tree["columns"]):
            label = tk.Label(self, text=col_name)
            label.grid(row=col_index, column=0, padx=5, pady=5)

            # entry_var = tk.StringVar()  çalışmadı, get boş dünüyor.
            entry = tk.Entry(self)
            entry.grid(row=col_index, column=1, padx=5, pady=5)
            entry.insert(0, selected_values[col_index])
            if col_name in self.disabled_colns:
                entry.config(state="readonly")
            self.entry_widgets.append(entry)  # Store entry widgets in the new list

        # Add a "Save" button
        save_button = tk.Button(self, text="Save", command=self.save_changes)
        save_button.grid(row=len(tree["columns"]), columnspan=2, pady=10)

    def disable_menu(self):
        pass  # no menu exist

    def enable_menu(self):
        pass  # no menu exist

    # Function to save the new values and update the Treeview
    def save_changes(self):
        new_values = [entry.get() for entry in self.entry_widgets]
        if len(new_values[1]) != 1:
            TKMessageBox("Error", "Letter must be length of one character !", messagebox.showerror, self).show()
        else:
            self.tree.item(self.selected_item, values=new_values)
            print_all_table_items(self.tree)
            self.destroy()


class AddNewRow(ModalWindow):
    def __init__(self, parent, tree):
        self.parent = parent
        self.tree = tree
        super().__init__(self.parent, resizable=(True, False), min_size=(200, 100))

        # Create a new window for inserting
        self.title("Insert Row")

        # new row index
        selected_item = self.tree.selection()
        if selected_item:
            selected_index = self.tree.index(selected_item[0])
            self.insert_pos = selected_index + 1
        else:
            self.insert_pos = "end"

        # Configure grid layout to allow resizing
        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)

        # Add labels and entry widgets for each column
        self.entry_widgets = []
        for col_index, col_name in enumerate(self.tree["columns"]):
            label = tk.Label(self, text=col_name)
            label.grid(row=col_index, column=0, padx=5, pady=5, sticky="w")

            entry = tk.Entry(self)
            entry.grid(row=col_index, column=1, padx=5, pady=5, sticky="we")
            self.entry_widgets.append(entry)

        # Add a "Save" button
        save_button = tk.Button(self, text="Save", command=self.save_changes)
        save_button.grid(row=len(self.tree["columns"]), columnspan=2, pady=10)

    # Function to save the new values and update the Treeview
    def save_changes(self):
        new_values = [entry.get() for entry in self.entry_widgets]
        self.tree.insert("", self.insert_pos, values=new_values)
        print_all_table_items(self.tree)
        self.destroy()


class AddNewRowLetterAnnot(ModalWindow):
    def __init__(self, parent, tree):
        self.parent = parent
        self.tree = tree
        super().__init__(self.parent, resizable=False)

        # Create a new window for inserting
        self.title("Insert Row")

        # new row index
        selected_item = self.tree.selection()
        if selected_item:
            selected_index = self.tree.index(selected_item[0])
            new_index = selected_index + 2  # +2 because indices start at 0
            self.insert_pos = selected_index + 1
        else:
            new_index = len(tree.get_children()) + 1
            self.insert_pos = "end"

        # Add labels and entry widgets for each column
        self.entry_widgets = []
        for col_index, col_name in enumerate(self.tree["columns"]):
            label = tk.Label(self, text=col_name)
            label.grid(row=col_index, column=0, padx=5, pady=5)

            entry = tk.Entry(self)
            entry.grid(row=col_index, column=1, padx=5, pady=5)
            if col_name == "Index":
                entry.insert(0, new_index)
                entry.config(state="readonly")
            self.entry_widgets.append(entry)

        # Add a "Save" button
        save_button = tk.Button(self, text="Save", command=self.save_changes)
        save_button.grid(row=len(self.tree["columns"]), columnspan=2, pady=10)

    # Function to save the new values and update the Treeview
    def save_changes(self):
        new_values = [entry.get() for entry in self.entry_widgets]
        if len(new_values[1]) != 1:
            TKMessageBox("Error", "Letter must be length of one character !", messagebox.showerror, self).show()
        else:
            new_item = self.tree.insert("", self.insert_pos, values=new_values)
            self.tree.selection_set(new_item)  # Select the newly added row
            update_indices(self.tree)
            print_all_table_items(self.tree)
            self.destroy()

    def enable_menu(self):
        pass

    def disable_menu(self):
        pass


def delete_selected_row(tree):
    selected_item = tree.selection()[0]  # Get the ID of the selected item
    if selected_item:
        tree.delete(selected_item)


def update_indices(tree):  # updates index column starting from 1
    for idx, item in enumerate(tree.get_children(), start=1):
        tree.item(item, values=(idx, *tree.item(item, 'values')[1:]))


class AddNewColumn(ModalWindow):
    def __init__(self, parent, tree):
        self.parent = parent
        self.tree = tree
        super().__init__(self.parent, resizable=(True, False), min_size=(200, 100))

        self.title("Add Column Name")
        self.geometry("200x100")

        label = tk.Label(self, text="Enter Column Name:")
        label.pack()

        self.entry = tk.Entry(self)
        self.entry.pack(padx=5, fill="x")

        add_button = tk.Button(self, text="Add", command=self.add_column)
        add_button.pack()

        self.mainloop()

    def add_column(self):
        # get coln name
        new_column_name = self.entry.get()

        # Preserve current column headers and their settings
        new_columns = [new_column_name]
        current_columns = list(self.tree['columns'])
        current_columns = {key: self.tree.heading(key) for key in current_columns}  # get heading info

        if new_column_name in current_columns:
            TKMessageBox("Error", "Column name aldready exist !", messagebox.showerror, self.parent).show()
            return

        # Update with new columns
        self.tree['columns'] = list(current_columns.keys()) + list(new_columns)
        for key in new_columns:
            self.tree.heading(key, text=key)

        # Set saved column values for the already existing columns
        for key in current_columns:
            # State is not valid to set with heading
            state = current_columns[key].pop('state')
            self.tree.heading(key, **current_columns[key])

        # Set default values for the new column for all existing items
        for item_id in self.tree.get_children():
            self.tree.set(item_id, new_column_name, "")

        print_all_table_items(self.tree)

        # destroy opened window
        self.destroy()


class DeleteColumn(ModalWindow):
    def __init__(self, parent, tree, current_columns):
        self.parent = parent
        self.tree = tree
        self.current_columns = current_columns
        super().__init__(self.parent, resizable=False)

        self.title("Delete Column")
        self.geometry("200x100")

        label = tk.Label(self, text="Column Name:")
        label.pack()

        self.col_name_var = tk.StringVar(self)
        self.col_name_var.set(self.current_columns[0])
        self.col_name_menu = tk.OptionMenu(self, self.col_name_var, *self.current_columns)
        self.col_name_menu.pack()

        add_button = tk.Button(self, text="Delete", command=self.delete_column)
        add_button.pack()

        self.mainloop()

    def delete_column(self):
        col_name = self.col_name_var.get()
        cids = self.tree["columns"]
        colns = [self.tree.heading(cid)["text"] for cid in cids]
        column_index = colns.index(col_name)

        # get heading info
        current_headings = {key: self.tree.heading(key) for index, key in enumerate(cids) if index != column_index}
        self.tree["columns"] = list(current_headings.keys())

        # Set saved column values for the already existing columns
        for key in current_headings:
            # State is not valid to set with heading
            _ = current_headings[key].pop('state')
            self.tree.heading(key, **current_headings[key])

        # Set default values for the new column for all existing items
        for item_id in self.tree.get_children():
            current_values = self.tree.item(item_id, 'values')
            updated_values = current_values[:column_index] + current_values[column_index + 1:]
            self.tree.item(item_id, values=updated_values)

        print_all_table_items(self.tree)

        # destroy opened window
        self.destroy()


def position_to_str(position, is_start):
    """Convert a Position object to its string representation."""
    diff = 1 if is_start else 0
    if isinstance(position, (ExactPosition, BeforePosition, AfterPosition, UncertainPosition)):
        return str(int(position)+diff)
    elif isinstance(position, WithinPosition) or isinstance(position, BetweenPosition):
        return f"{position._left+diff},{position._right+diff}"
    elif isinstance(position, OneOfPosition):
        return ",".join(map(lambda x: str(x+diff), position.position_choices))
    else:
        return ""


def loc_str(location):
    answer = f"[{position_to_str(location._start, is_start=True)}:{position_to_str(location._end, is_start=False)}]"
    if location.ref and location.ref_db:
        answer = f"{location.ref_db}:{location.ref}{answer}"
    elif location.ref:
        answer = location.ref + answer
    # Is ref_db without ref meaningful?
    if location.strand is None:
        return answer
    elif location.strand == +1:
        return answer + "(+)"
    elif location.strand == -1:
        return answer + "(-)"
    else:
        # strand = 0, stranded but strand unknown, ? in GFF3
        return answer + "(?)"


class NestedObjectViewer:
    def __init__(self, parent_window, tab, data):
        self.parent_window = parent_window
        self.data = data

        # Treeview with both tree and headings
        self.tree = ttk.Treeview(tab, columns=("Value",), show="tree headings")
        self.tree.heading("#0", text="Field")  # Tree column (for keys)
        self.tree.heading("Value", text="Value")  # Value column
        self.tree.column("#0", width=300, anchor="w")  # Tree column
        self.tree.column("Value", width=400, anchor="w")  # Value column
        # self.tree.pack(fill=tk.BOTH, expand=True)

        # Create vertical scrollbar
        vscrollbar = ttk.Scrollbar(tab, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vscrollbar.set)
        vscrollbar.pack(side="right", fill="y")

        # Create horizontal scrollbar
        hscrollbar = ttk.Scrollbar(tab, orient="horizontal", command=self.tree.xview)
        self.tree.configure(xscrollcommand=hscrollbar.set)
        hscrollbar.pack(side="bottom", fill="x")

        # Bind the Configure event to adjust column widths
        # self.tree.bind("<Configure>", lambda event: self.tree_on_configure(event, self.tree))

        # Populate the Treeview
        self.data_map = {}  # To map row IDs to object fields
        self.populate_tree(self.data)

        # Bind selection event
        self.tree.bind("<<TreeviewSelect>>", self.on_row_selected)

    def populate_tree(self, data, parent="", data_info = None):  # data_info: data_insert_index, container, key
        """Recursively populate the Treeview with data."""

        # add row for data  (if enabled)
        if data_info:
            data_insert_index, container, key = data_info
            value = "" if isinstance(data, (dict, list, Reference, SeqFeature)) else str(data)
            node_id = self.tree.insert(parent, data_insert_index, text=key, values=(value,))
            self.data_map[node_id] = (container, (key, data))  # Store in data_map
            parent = node_id  # insert children under data row

        # add rows for children of data
        if isinstance(data, dict):
            self._process_iterable(data.items(), parent, container=data, key_func=lambda k_v: k_v[0], value_func=lambda k_v: k_v[1])
        elif isinstance(data, list):
            self._process_iterable(enumerate(data), parent, container=data, key_func=lambda i_v: f"[{i_v[0]+1}]", value_func=lambda i_v: i_v[1])
        elif isinstance(data, (Reference, SeqFeature)):
            self._process_iterable(vars(data).items(), parent, container=data, key_func=lambda k_v: k_v[0], value_func=lambda k_v: k_v[1])

    def _process_iterable(self, iterable, parent, container, key_func, value_func):
        """Helper method to process dictionaries, lists, and custom objects."""
        for item in iterable:
            key = key_func(item)
            value = value_func(item)

            if isinstance(value, (dict, list, Reference, SeqFeature)):
                # Insert parent node for nested structures or custom objects
                node_id = self.tree.insert(parent, "end", text=key, values=("",))
                self.data_map[node_id] = (container, item)  # Store in data_map
                self.populate_tree(value, node_id)  # Recursively populate
            else:
                # Primitive value or Location: key and value in the same row
                str_fun = loc_str if isinstance(value, Location) else str
                node_id = self.tree.insert(parent, "end", text=key, values=(str_fun(value),))
                self.data_map[node_id] = (container, item)

    def on_row_selected(self, event):
        """Handle row selection and find the corresponding field."""
        if len(self.tree.selection()) > 0:
            selected_item = self.tree.selection()[0]  # Get selected row ID
            field_info = self.data_map[selected_item]  # Retrieve field info
            print("Selected:", field_info)

    def edit_selected_row(self):
        selected_item = self.tree.selection()[0]  # Get selected row ID
        container, item = self.data_map[selected_item]  # Retrieve field info

        # Location cases
        if isinstance(item[1], Location):  # edit simple loc
            self.open_location_window_edit(selected_item, item, container)
        elif item[1] is None and item[0] == "location" and isinstance(container, SeqFeature):  # editing loc after adding a new SeqFeature
            s = self.LocationTypeSelector(self.parent_window)
            if s.is_closed_succesfully:
                if s.selected_option.get() == "simple":
                    self.open_location_window_edit(selected_item, item, container)
                else:  # compound
                    # update table
                    self.tree.item(selected_item, text="location", values=("",))

                    # update data-table map
                    self.data_map[selected_item] = (container, (item[0], []))  # list item

        else:  # primitive
            is_field_disabled = not isinstance(container, dict)  # disabled if not dict item (e.g list item, object field)
            is_value_disabled = isinstance(item[1], (Reference, SeqFeature, dict, list))  # enabled if primitive

            if is_field_disabled and is_value_disabled:
                return  # nothing to edit

            field = f"[{item[0] + 1}]" if isinstance(container, list) else item[0]
            value = item[1] if not is_value_disabled else ""

            fields_reserved = [key for key in container if key != field] if isinstance(container, dict) else None

            w = self.FieldValueWindow(self.parent_window, self.tree, window_title="Edit Selected Row", field=field, value=value,
                                      is_field_disabled=is_field_disabled, is_value_disabled=is_value_disabled, fields_reserved=fields_reserved)
            if w.is_closed_succesfully:
                print(f"Saved Field: {w.field}, Value: {w.value}")

                # update data
                key = item[0] if isinstance(container, list) else w.field
                value = w.value if not is_value_disabled else item[1]  # old value if disabled
                if isinstance(container, list):
                    container[key] = value  # cannot use w.field as it is like "[item[0]]"
                elif isinstance(container, dict):
                    del container[field]
                    container[key] = value
                else:  # object
                    setattr(container, key, value)

                # update table
                self.tree.item(selected_item, text=w.field, values=(w.value,))

                # update data-table map
                self.data_map[selected_item] = (container, (key, value))  # (container,item)

    def open_location_window_edit(self, selected_item, item, container):
        w = self.SimpleLocationWindow(self.parent_window, item[1])
        if w.is_closed_succesfully:
            # update data
            key = item[0]
            value = w.location
            if isinstance(container,
                          list):  # parent is list for "compound loc" or list for "location" attribute of Reference object
                container[key] = value
            else:  # object (SeqFeature)
                setattr(container, key, value)

            # update table
            self.tree.item(selected_item, text=f"[{item[0] + 1}]" if isinstance(container, list) else key,
                           values=(loc_str(value),))

            # update data-table map
            self.data_map[selected_item] = (container, (key, value))  # (container,item)

    def delete_selected_row(self):
        selected_item = self.tree.selection()[0]  # Get selected row ID
        container, item = self.data_map[selected_item]  # Retrieve field info

        if not isinstance(container, (list, dict)):  # only delete list and dict items
            return

        # Delete the selected item from table
        parent = self.tree.parent(selected_item)  # Get parent row ID
        self.tree.delete(selected_item)

        # Delete the selected item from data
        del container[item[0]]  # key or index

        # update data-table map
        del self.data_map[selected_item]

        # reorder if list item deleted
        if isinstance(container, list):
            self.reoorder_list_items(parent)

    def add_new_row(self, is_annot_tab):

        if len(self.tree.selection()) == 0:  # add top level row
            if is_annot_tab:  # add str or list of str (into annot dict)
                w = self.RowTypeSelector(self.parent_window)
                value_type = str if w.selected_option.get() == "single" else list
                if not w.is_closed_succesfully:
                    return

                container = self.data  # record.annotations
                children = list(self.tree.get_children(""))
                fields_reserved = [key for key in container]

                w = self.FieldValueWindow(self.parent_window, self.tree, window_title="Add New Row", field="",
                                          value="",
                                          is_field_disabled=False, is_value_disabled=(value_type == list),
                                          fields_reserved=fields_reserved)
                if not w.is_closed_succesfully:
                    return

                # create key and value
                key = w.field
                value = list() if value_type == list else w.value

                # add to data
                container[key] = value

                # add into table
                data_info = len(children), container, key
                self.populate_tree(value, parent="", data_info=data_info)

            else:  # add Feature row (into Feature list)
                container = self.data  # record.features
                children = list(self.tree.get_children(""))

                # create key and value
                key = len(children)
                value = SeqFeature()

                # add to data
                container.insert(key, value)

                # add into table
                data_info = key, container, key
                self.populate_tree(value, parent="", data_info=data_info)

                # reorder if list item
                self.reoorder_list_items("")

        else:  # row will be added after selected one
            selected_item = self.tree.selection()[0]  # Get selected row ID
            container, sel_item = self.data_map[selected_item]  # Retrieve field info

            # Location cases
            # check if parent is "location" of Reference object
            parent_id = self.tree.parent(selected_item)
            if parent_id != "":
                parent_container, parent_item = self.data_map[parent_id]  # Retrieve field info
                if parent_item[0] == "location" and isinstance(parent_container, Reference):
                    self.add_new_row_location(selected_item, parent_id, parent_container, parent_item)
                    return

            # check if selected is Location object, parent is list(compound) and we are on "Annot" tab
            if isinstance(sel_item[1], Location) and isinstance(container, list) and is_annot_tab:
                parent_container, parent_item = self.data_map[parent_id]  # Retrieve field info
                if parent_item[0] != "location": # make sure it is not "location" list in Reference, but a compound loc list
                    children = list(self.tree.get_children(parent_id))  # Convert to list for indexing
                    insert_position = children.index(selected_item) + 1  # Insert after sel row
                    w = self.SimpleLocationWindow(self.parent_window, None)
                    if w.is_closed_succesfully:
                        self.add_simple_into_compound(insert_position, w.location, container, parent_id)
                    return

            # check if selected is Location object, parent is list(compound) and we are on "Features" tab
            if isinstance(sel_item[1], Location) and isinstance(container, list) and not is_annot_tab:
                children = list(self.tree.get_children(parent_id))  # Convert to list for indexing
                insert_position = children.index(selected_item) + 1  # Insert after sel row
                w = self.SimpleLocationWindow(self.parent_window, None)
                if w.is_closed_succesfully:
                    self.add_simple_into_compound(insert_position, w.location, container, parent_id)
                return

            if not isinstance(container, (list, dict)):  # can only add into list and dict (not into object e.g.Feature)
                return

            is_field_disabled = not isinstance(container, dict)  # disabled if not dict item
            is_value_disabled = isinstance(sel_item[1], (Reference, SeqFeature, dict, list))  # enabled if primitive

            new_field = ""
            new_value = ""
            if not (is_field_disabled and is_value_disabled):  # Open window if either is editable
                field = f"[{sel_item[0] + 2}]" if isinstance(container, list) else ""
                value = ""

                fields_reserved = [key for key in container] if isinstance(container, dict) else None

                w = self.FieldValueWindow(self.parent_window, self.tree, window_title="Add New Row", field=field,
                                          value=value,
                                          is_field_disabled=is_field_disabled, is_value_disabled=is_value_disabled,
                                          fields_reserved=fields_reserved)
                if w.is_closed_succesfully:
                    print(f"Saved Field: {w.field}, Value: {w.value}")
                    new_field = w.field
                    new_value = w.value
                else:
                    return

            # create key and value
            key = sel_item[0]+1 if isinstance(container, list) else new_field
            if isinstance(sel_item[1], Reference):
                value = Reference()
            elif isinstance(sel_item[1], SeqFeature):
                value = SeqFeature()
            elif isinstance(sel_item[1], list):
                value = list()
            elif isinstance(sel_item[1], dict):
                value = dict()
            else:
                value = new_value  # primitive

            # add to data
            if isinstance(container, list):
                container.insert(key, value)
            elif isinstance(container, dict):
                container[key] = value
            else:  # object
                setattr(container, key, value)

            # add into table
            parent = self.tree.parent(selected_item)  # Get parent row ID
            children = list(self.tree.get_children(parent))  # Convert to list for indexing
            selected_index = children.index(selected_item)
            insert_position = selected_index + 1  # Insert after sel row
            data_info = insert_position, container, key
            self.populate_tree(value, parent=parent, data_info=data_info)

            # reorder if list item
            if isinstance(container, list):
                self.reoorder_list_items(parent)

    def add_new_row_location(self, selected_item, parent_id, parent_container, parent_item):
        container = parent_item[1]
        children = list(self.tree.get_children(parent_id))  # Convert to list for indexing

        s = self.LocationTypeSelector(self.parent_window)
        if s.is_closed_succesfully:
            if s.selected_option.get() == "simple":
                w = self.SimpleLocationWindow(self.parent_window, None)
                if w.is_closed_succesfully:
                    # update data
                    key = children.index(selected_item) + 1  # Insert after sel row
                    value = w.location

                    # add to data
                    container.insert(key, value)

                    # add into table
                    node_id = self.tree.insert(parent_id, key, text=f"[{key + 1}]", values=(loc_str(value),))

                    # update data-table map
                    self.data_map[node_id] = (container, (key, value))  # (container,item)

                    # reorder if list item
                    self.reoorder_list_items(parent_id)
            else:  # compound
                key = children.index(selected_item) + 1  # Insert after sel row
                value = []

                # add to data
                container.insert(key, value)

                # add into table
                node_id = self.tree.insert(parent_id, key, text=f"[{key + 1}]", values=("",))

                # update data-table map
                self.data_map[node_id] = (container, (key, value))  # (container,item)

                # reorder if list item
                self.reoorder_list_items(parent_id)

    def add_sub_row(self, is_annot_tab):  # add into first index

        if len(self.tree.selection()) == 0:
            return

        # new row will be added under the selected one (container is selected one)
        selected_item = self.tree.selection()[0]  # Get selected row ID
        sel_parent, (sel_key, container) = self.data_map[selected_item]  # Retrieve field info

        # Location cases
        if sel_key == "location" and isinstance(sel_parent, Reference): # add (simple or compound) into location field of Reference
            self.add_sub_row_location(selected_item, sel_parent, sel_key, container)
            return
        elif isinstance(container, list):  # check if selected is "compound loc" in Reference object location list.
            parent_id = self.tree.parent(selected_item)
            if parent_id != "":
                parent_container, parent_item = self.data_map[parent_id]  # Retrieve field info
                if parent_item[0] == "location" and isinstance(parent_container, Reference):
                    w = self.SimpleLocationWindow(self.parent_window, None)
                    if w.is_closed_succesfully:
                        self.add_simple_into_compound(0, w.location, container, selected_item)
                    return
        elif sel_key == "location" and isinstance(sel_parent, SeqFeature) and isinstance(container, list): # SeqFeature location is compound
            w = self.SimpleLocationWindow(self.parent_window, None)
            if w.is_closed_succesfully:
                self.add_simple_into_compound(0, w.location, container, selected_item)
            return

        # define "value type" for item to be inserted
        children = list(self.tree.get_children(selected_item))
        item_value_type = str
        if len(children) > 0 :
            (_, (_, first_child)) = self.data_map[children[0]]
            item_value_type = type(first_child)  # same as child type
        elif isinstance(sel_parent, SeqFeature) and sel_key == "qualifiers":  # Features tab
            item_value_type = list
        elif sel_key == "references" and is_annot_tab:   # Annotations tab
            item_value_type = Reference

        # can only add into list and dict (not into object e.g.Feature)
        if not isinstance(container, (list, dict)):
            return

        is_field_disabled = not isinstance(container, dict)  # disabled if not dict item
        is_value_disabled = any(item_value_type == t for t in (Reference, SeqFeature, dict, list)) # enabled if primitive

        new_field = ""
        new_value = ""
        if not (is_field_disabled and is_value_disabled):  # Open window if either is editable
            field = f"[1]" if isinstance(container, list) else ""  # insert at top
            value = ""

            fields_reserved = [key for key in container] if isinstance(container, dict) else None

            w = self.FieldValueWindow(self.parent_window, self.tree, window_title="Add Sub-Row", field=field,
                                      value=value,
                                      is_field_disabled=is_field_disabled, is_value_disabled=is_value_disabled,
                                      fields_reserved=fields_reserved)
            if w.is_closed_succesfully:
                print(f"Saved Field: {w.field}, Value: {w.value}")
                new_field = w.field
                new_value = w.value
            else:
                return

        # create key and value
        key = 0 if isinstance(container, list) else new_field
        if item_value_type == Reference:
            value = Reference()
        elif item_value_type == SeqFeature:
            value = SeqFeature()
        elif item_value_type == list:
            value = list()
        elif item_value_type == dict:
            value = dict()
        else:
            value = new_value  # primitive

        # add to data
        if isinstance(container, list):
            container.insert(key, value)
        elif isinstance(container, dict):
            container[key] = value
        else:  # object
            setattr(container, key, value)

        # add into table
        data_info = 0, container, key
        self.populate_tree(value, parent=selected_item, data_info=data_info)

        # reorder if list item
        if isinstance(container, list):
            self.reoorder_list_items(selected_item)

    def add_sub_row_location(self, selected_item, sel_parent, sel_key, container):  #  Reference -> location -> add simple/compound
        s = self.LocationTypeSelector(self.parent_window)
        if s.is_closed_succesfully:
            if s.selected_option.get() == "simple":
                w = self.SimpleLocationWindow(self.parent_window, None)
                if w.is_closed_succesfully:
                    # update data
                    key = 0
                    value = w.location

                    # add to data
                    container.insert(key, value)

                    # add into table
                    node_id = self.tree.insert(selected_item, 0, text="[1]", values=(loc_str(value),))

                    # update data-table map
                    self.data_map[node_id] = (container, (key, value))  # (container,item)

                    # reorder if list item
                    self.reoorder_list_items(selected_item)

            else:  # compound
                key = 0
                value = []

                # add to data
                container.insert(key, value)

                # add into table
                node_id = self.tree.insert(selected_item, key, text=f"[{key + 1}]", values=("",))

                # update data-table map
                self.data_map[node_id] = (container, (key, value))  # (container,item)

                # reorder if list item
                self.reoorder_list_items(selected_item)

    def add_simple_into_compound(self, index, location, container, parent_id):
        # update data
        key = index
        value = location

        # add to data
        container.insert(key, value)

        # add into table
        node_id = self.tree.insert(parent_id, index, text=f"[{index+1}]", values=(loc_str(value),))

        # update data-table map
        self.data_map[node_id] = (container, (key, value))  # (container,item)

        # reorder if list item (required if add into top)
        self.reoorder_list_items(parent_id)

    def add_special_row(self, field):
        container = self.data  # record.annotations
        children = list(self.tree.get_children(""))
        fields_reserved = [key for key in container]
        if field in fields_reserved:
            TKMessageBox("Error", f"Field {field} is already used !", messagebox.showerror, self.parent_window).show()
            return

        value_type = special_field_types[field]
        new_value = ""
        if not value_type == list:
            w = self.FieldValueWindow(self.parent_window, self.tree, window_title="Add Genbank/EMBL Row", field=field,
                                      value=0 if value_type == int else "",
                                      is_field_disabled=True, is_value_disabled=False)
            if not w.is_closed_succesfully:
                return
            else:
                new_value = w.value

        # create key and value
        key = field
        value = list() if value_type == list else new_value

        # add to data
        container[key] = value

        # add into table
        data_info = len(children), container, key
        self.populate_tree(value, parent="", data_info=data_info)

    def reoorder_list_items(self, parent):
        # Get remaining children under the same parent
        children = self.tree.get_children(parent)

        # Reorder text and update keys in data_map
        for new_index, child in enumerate(children):
            container, (_, value) = self.data_map[child]  # Extract container and value

            # Update the tree row text to "[index]" (starting from 1)
            current_values = self.tree.item(child, "values")  # Preserve value display
            self.tree.item(child, text=f"[{new_index + 1}]", values=current_values)

            # Update the index in data_map
            self.data_map[child] = (container, (new_index, value))

    class FieldValueWindow(ModalWindow):
        def __init__(self, parent, tree, window_title="", field="", value="", is_field_disabled=False, is_value_disabled=False, fields_reserved=None):
            self.parent = parent
            self.tree = tree
            self.columns = ["Field", "Value"]
            super().__init__(self.parent, min_size=(400, 200))

            self.title(window_title)
            self.geometry("400x200")
            self.is_closed_succesfully = False

            # Field label + Entry
            ttk.Label(self, text="Field:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
            self.field_entry = ttk.Entry(self)
            self.field_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=5)
            self.field_entry.insert(0, field)
            if is_field_disabled:  # Disable the entry widget if needed
                self.field_entry.config(state="readonly")
            self.fields_reserved = fields_reserved

            # Value label + Text
            ttk.Label(self, text="Value:").grid(row=1, column=0, sticky="nw", padx=5, pady=5)
            self.value_text = tk.Text(self, height=5, wrap="word")
            self.value_text.grid(row=1, column=1, sticky="nsew", padx=5, pady=5)
            self.value_text.insert("1.0", value if value is not None else "")
            if is_value_disabled: # Disable the text widget if needed
                self.value_text.config(state="disabled")
            self.value_type = type(value)

            # Save button centered
            self.save_button = ttk.Button(self, text="Save", command=self.on_save)
            self.save_button.grid(row=2, column=0, columnspan=2, pady=10)

            # Grid weight for resizing
            self.columnconfigure(1, weight=1)  # Allow horizontal resizing
            self.rowconfigure(1, weight=1)  # Allow vertical resizing for Text

            # Pause execution in the annot window until this window is closed
            self.parent.wait_window(self)

        def on_save(self):
            self.field = self.field_entry.get()
            if self.fields_reserved is not None and self.field in self.fields_reserved:
                TKMessageBox("Error", f"Field {self.field} is already used !", messagebox.showerror, self).show()
                return

            value = self.value_text.get("1.0", "end").strip()
            self.value = self.convert_to_original_type(value)

            self.is_closed_succesfully = True
            self.destroy()  # Close the window

        def convert_to_original_type(self, value_str):
            """Convert the string input back to its original type."""
            if self.value_type == int:
                try:
                    return int(value_str)
                except ValueError:
                    return value_str  # Keep as string if conversion fails
            elif self.value_type == float:
                try:
                    return float(value_str)
                except ValueError:
                    return value_str # Keep as string if conversion fails
            else:
                return value_str

        def disable_menu(self):
            pass

        def enable_menu(self):
            pass

    class RowTypeSelector(ModalWindow):
        def __init__(self, parent):
            self.parent = parent
            super().__init__(self.parent, resizable=False)
            self.title("Select Row Type")
            # self.geometry("300x150")

            self.is_closed_succesfully = False
            self.selected_option = tk.StringVar(value="single")  # Default: single row

            # Radio Buttons
            ttk.Radiobutton(self, text="Single Row", variable=self.selected_option, value="single").pack(
                anchor="w", padx=20, pady=10)
            ttk.Radiobutton(self, text="Row with Multiple Item", variable=self.selected_option,
                            value="list").pack(anchor="w", padx=20)

            # OK Button
            ttk.Button(self, text="OK", command=self.on_ok).pack(pady=10)

            # Pause execution until this window is closed
            self.parent.wait_window(self)

        def on_ok(self):
            """Return selected option and close the window."""
            print("Selected:", self.selected_option.get())  # Replace with return logic
            self.is_closed_succesfully = True
            self.destroy()

    class SimpleLocationWindow(ModalWindow):
        def __init__(self, parent, location=None):
            self.parent = parent
            super().__init__(self.parent, resizable=False)
            self.title("Location")
            pad = 5
            self.is_closed_succesfully = False
            self.location = location

            # Position Type Selection
            ttk.Label(self, text="Start Position:").pack(pady=pad, padx=pad, anchor="w")
            value = POSITION_VALUES[type(self.location.start).__name__] if self.location is not None else POSITION_VALUES["ExactPosition"]
            self.start_pos_type = tk.StringVar(value=value)
            pos_options = list(POSITION_VALUES.values())
            self.start_pos_menu = tk.OptionMenu(self, self.start_pos_type, *pos_options, command=self.update_start_hint)
            self.start_pos_menu.pack(fill="x", padx=10)

            # Start Position Entry with Placeholder
            self.start_entry = ttk.Entry(self, foreground="gray")
            self.start_hint = POSITION_HINTS["ExactPosition"]
            self.start_entry.insert(0, self.start_hint)
            if self.location is not None:
                self.start_entry.delete(0, 'end')
                self.start_entry.config(foreground="black")
                self.start_entry.insert(0, position_to_str(self.location.start, True))
            self.start_entry.pack(fill="x", padx=10, pady=2)

            # Bind events
            self.start_entry.bind("<Button-1>", self.click_start)
            self.start_entry.bind("<Leave>", self.leave_start)

            # End Position Type
            ttk.Label(self, text="End Position:").pack(pady=pad, padx=pad, anchor="w")
            value = POSITION_VALUES[type(self.location.end).__name__] if self.location is not None else POSITION_VALUES["ExactPosition"]
            self.end_pos_type = tk.StringVar(value=value)
            self.end_pos_menu = tk.OptionMenu(self, self.end_pos_type, *pos_options, command=self.update_end_hint)
            self.end_pos_menu.pack(fill="x", padx=10)

            # End Position Entry with Placeholder
            self.end_entry = ttk.Entry(self, foreground="gray")
            self.end_hint = POSITION_HINTS["ExactPosition"]
            self.end_entry.insert(0, self.end_hint)
            if self.location is not None:
                self.end_entry.delete(0, 'end')
                self.end_entry.config(foreground="black")
                self.end_entry.insert(0, position_to_str(self.location.end, False))
            self.end_entry.pack(fill="x", padx=10, pady=2)

            # Bind events
            self.end_entry.bind("<Button-1>", self.click_end)
            self.end_entry.bind("<Leave>", self.leave_end)

            # Strand Selection
            ttk.Label(self, text="Strand:").pack(pady=pad, padx=pad, anchor="w")
            strand_map = {1: "Forward", -1: "Reverse", 0: "Unknown", None: "None"}
            self.strand = tk.StringVar(value=strand_map[self.location.strand] if self.location is not None else "None")
            strand_options = ["Forward", "Reverse", "Unknown", "None"]
            self.strand_menu = tk.OptionMenu(self, self.strand, *strand_options)
            self.strand_menu.pack(fill="x", padx=10)

            # Reference & Database Fields
            ttk.Label(self, text="Reference:").pack(pady=pad, padx=pad, anchor="w")
            self.ref_entry = ttk.Entry(self)
            self.ref_entry.insert(0, location.ref if (location is not None and location.ref is not None) else "")
            self.ref_entry.pack(fill="x", padx=10)

            ttk.Label(self, text="Reference Database:").pack(pady=pad, padx=pad, anchor="w")
            self.ref_db_entry = ttk.Entry(self)
            self.ref_db_entry.insert(0, location.ref_db if (location is not None and location.ref_db is not None) else "")
            self.ref_db_entry.pack(fill="x", padx=10)

            # OK Button
            ttk.Button(self, text="OK", command=self.on_ok).pack(pady=10)

            # Pause execution in the annot window until this window is closed
            self.parent.wait_window(self)

        # Click event for start position entry
        def click_start(self, *args):
            if self.start_entry.get().strip() == self.start_hint:
                self.start_entry.delete(0, 'end')
                self.start_entry.config(foreground="black")

        # Leave event for start position entry
        def leave_start(self, *args):
            if self.start_entry.get().strip() == "":
                self.start_entry.insert(0, self.start_hint)
                self.start_entry.config(foreground="gray")
            self.focus()

        # Click event for end position entry
        def click_end(self, *args):
            if self.end_entry.get().strip() == self.end_hint:
                self.end_entry.delete(0, 'end')
                self.end_entry.config(foreground="black")

        # Leave event for end position entry
        def leave_end(self, *args):
            if self.end_entry.get().strip() == "":
                self.end_entry.insert(0, self.end_hint)
                self.end_entry.config(foreground="gray")
            self.focus()

        # Update hints dynamically when type changes
        def update_start_hint(self, *args):
            self.start_hint = POSITION_HINTS[self.start_pos_type.get().replace(" ","")]
            self.start_entry.delete(0, tk.END)
            self.start_entry.insert(0, self.start_hint)
            self.start_entry.config(foreground="gray")
            if self.start_pos_type.get() == "UnknownPosition":
                self.start_entry.config(state="disabled")
            else:
                self.start_entry.config(state="normal")
                self.start_hint = POSITION_HINTS[self.start_pos_type.get().replace(" ","")]
                self.start_entry.delete(0, tk.END)
                self.start_entry.insert(0, self.start_hint)
                self.start_entry.config(foreground="gray")

        def update_end_hint(self, *args):
            self.end_hint = POSITION_HINTS[self.end_pos_type.get().replace(" ","")]
            self.end_entry.delete(0, tk.END)
            self.end_entry.insert(0, self.end_hint)
            self.end_entry.config(foreground="gray")
            if self.end_pos_type.get() == "UnknownPosition":
                self.end_entry.config(state="disabled")
            else:
                self.end_entry.config(state="normal")
                self.end_hint = POSITION_HINTS[self.end_pos_type.get().replace(" ","")]
                self.end_entry.delete(0, tk.END)
                self.end_entry.insert(0, self.end_hint)
                self.end_entry.config(foreground="gray")

        def parse_position(self, pos_type, pos_value, is_start):  # -1 for python indexing that Position objects require
            """Parse position input and create appropriate Position object."""
            try:

                if any(pos_type == pos for pos in ["Exact Position", "Before Position", "After Position", "Uncertain Position"]):
                    pos_value = int(pos_value) - 1 if is_start else int(pos_value)

                    if pos_type == "Exact Position":
                        return ExactPosition(pos_value)
                    elif pos_type == "Before Position":
                        return BeforePosition(pos_value)
                    elif pos_type == "After Position":
                        return AfterPosition(pos_value)
                    elif pos_type == "Uncertain Position":
                        return UncertainPosition(pos_value)

                if any(pos_type == pos for pos in ["Within Position", "Between Position"]):
                    left, right = map(int, pos_value.split(","))
                    assert left <= right
                    if is_start:
                        left -= 1
                        right -= 1
                        default = left
                    else:
                        default = right

                    if pos_type == "Within Position":
                        return WithinPosition(default, left, right)
                    elif pos_type == "Between Position":
                        return BetweenPosition(default, left, right)

                elif pos_type == "One Of Position":
                    choices = list(map(int, pos_value.split(",")))
                    diff = 1 if is_start else 0
                    choices = list(map(lambda x: x-diff, choices))
                    return OneOfPosition(choices[0], choices)

                elif pos_type == "Unknown Position":
                    return UnknownPosition()
            except Exception:
                pos_type = pos_type.replace(" ", "")
                TKMessageBox("Invalid Input", f"Invalid format for {pos_type}. Expected: {POSITION_HINTS[pos_type]}", messagebox.showerror, self.parent).show()
                return None

        def on_ok(self):
            """Create SimpleLocation when OK is pressed."""
            start_pos = self.parse_position(self.start_pos_type.get(), self.start_entry.get(), is_start=True)
            if start_pos is None :
                return  # Stop if invalid input

            end_pos = self.parse_position(self.end_pos_type.get(), self.end_entry.get(), is_start=False)
            if end_pos is None:
                return  # Stop if invalid input

            strand_map = {"Forward": 1, "Reverse": -1, "Unknown": 0, "None": None}
            strand_value = strand_map[self.strand.get()]

            ref = self.ref_entry.get()
            ref_db = self.ref_db_entry.get()

            try:
                self.location = SimpleLocation(start_pos, end_pos, strand=strand_value,
                                               ref=ref if ref != "" else None,
                                               ref_db=ref_db if ref_db != "" else None)
            except ValueError as e:
                TKMessageBox("Invalid Input", e,
                             messagebox.showerror, self.parent).show()
                return # Stop if invalid input

            print("Created SimpleLocation:", self.location)  # Replace with return logic
            self.is_closed_succesfully = True
            self.destroy()

    class LocationTypeSelector(ModalWindow):
        def __init__(self, parent):
            self.parent = parent
            super().__init__(self.parent, resizable=False)
            self.title("Select Location Type")
            # self.geometry("300x150")

            self.is_closed_succesfully = False
            self.selected_option = tk.StringVar(value="simple")  # Default: single row

            # Radio Buttons
            ttk.Radiobutton(self, text="Simple Location", variable=self.selected_option, value="simple").pack(
                anchor="w", padx=20, pady=10)
            ttk.Radiobutton(self, text="Compound Location", variable=self.selected_option,
                            value="compound").pack(anchor="w", padx=20)

            # OK Button
            ttk.Button(self, text="OK", command=self.on_ok).pack(pady=10)

            # Pause execution until this window is closed
            self.parent.wait_window(self)

        def on_ok(self):
            """Return selected option and close the window."""
            print("Selected:", self.selected_option.get())  # Replace with return logic
            self.is_closed_succesfully = True
            self.destroy()
