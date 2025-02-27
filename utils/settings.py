from utils.model_window import ModalWindow
import os
import json
import tkinter as tk
from tkinter import ttk, filedialog
from appdirs import AppDirs
from filelock import FileLock, Timeout
from common.configs import APP_NAME
from utils.utils import check_font_size_range
from utils.windows import TKMessageBox
import webbrowser

# Set up app-specific directories
dirs = AppDirs(APP_NAME)
os.makedirs(dirs.user_data_dir, exist_ok=True)  # Ensure directory exists
# (Modules are run only once even if they are imported multiple times)

# paths for settings
settings_general_path = os.path.join(dirs.user_data_dir, "settings_general.json")
settings_tools_path = os.path.join(dirs.user_data_dir, "settings_tools.json")
settings_othertools_path = os.path.join(dirs.user_data_dir, "settings_othertools.json")
print(settings_general_path, settings_tools_path, settings_othertools_path, sep="\n")

# Define the lock file paths for file lock
settings_general_lock_path = settings_general_path + ".lock"
settings_tools_lock_path = settings_tools_path + ".lock"
settings_othertools_lock_path = settings_othertools_path + ".lock"

# Define tools and their specific items for each tab
# For format names, we stick to the related website usage.
tool_params = {
    "Alignment Tools": {
        'Muscle': [('Executable Path', 'entry_browse'), ('Input File', 'entry_browse'),
                   ('Output File', 'entry_browse'), ('Output Format', 'optionmenu', ["FASTA"]),
                   ('Tool Website', 'button', 'https://drive5.com/muscle5/'),
                   ('Version Checked', 'label', 'v5.2'),
                   ],
        'ClustalOmega': [('Executable Path', 'entry_browse'), ('Input File', 'entry_browse'),
                         ('Output File', 'entry_browse'), ('Output Format', 'optionmenu',
                                                           ["fasta", "clustal", "msf", "phylip", "selex",
                                                            "stockholm", "vienna"]),
                         ('Tool Website', 'button', 'http://www.clustal.org/omega/'),
                         ('Version Checked', 'label', 'v1.2.3')],
        'Prank': [('Executable Path', 'entry_browse'), ('Input File', 'entry_browse'),
                  ('Output File', 'entry_browse'), ('Output Format', 'optionmenu',
                                                    ["IG/Stanford", "Pearson/Fasta", "GenBank/GB", "Phylip3.2",
                                                     "NBRF", "Phylip", "EMBL", "PIR/CODATA", "DNAStrider",
                                                     "MSF", "Fitch", "PAUP/NEXUS"]),
                  ('Tool Website', 'button', 'https://ariloytynoja.github.io/prank-msa/'),
                  ('Version Checked', 'label', 'v170427')],
        'Mafft': [('Executable Path', 'entry_browse'), ('Input File', 'entry_browse'),
                  ('Output File', 'entry_browse'),
                  ('Output Format', 'optionmenu', ["fasta", "phylip", "clustal"]),
                  ('Tool Website', 'button', 'https://mafft.cbrc.jp/alignment/software/'),
                  ('Version Checked', 'label', 'version 7')],
        'Probcons': [('Executable Path', 'entry_browse'), ('Input File', 'entry_browse'),
                     ('Output File', 'entry_browse'), ('Output Format', 'optionmenu', ["MFA", "ClustalW"]),
                     ('Tool Website', 'button', 'http://probcons.stanford.edu/'),
                     ('Version Checked', 'label', 'version 1.12')],
        'MSAProbs': [('Executable Path', 'entry_browse'), ('Input File', 'entry_browse'),
                     ('Output File', 'entry_browse'), ('Output Format', 'optionmenu', ["FASTA", "CLUSTALW"]),
                     ('Tool Website', 'button', 'https://msaprobs.sourceforge.net/homepage.htm#latest'),
                     ('Version Checked', 'label', 'v0.9.7')],
        'Needleman-Wunsch (Global)': [('Executable Path', 'entry_browse'), ('First Input File', 'entry_browse'),
                                      ('Second Input File', 'entry_browse'), ('Gap Open', 'entry'),
                                      ('Gap Extend', 'entry'), ('Output File', 'entry_browse'), (
                                          'Output Format', 'optionmenu',
                                          ["markx0", "markx1", "markx2", "markx3", "markx10", "srspair", "score"]),
                                      ('Tool Website', 'button', 'https://emboss.sourceforge.net/download/'),
                                      ('Version Checked', 'label', 'EMBOSS-6.6.0')],
        'Smith-Waterman (Local)': [('Executable Path', 'entry_browse'), ('First Input File', 'entry_browse'),
                                   ('Second Input File', 'entry_browse'), ('Gap Open', 'entry'),
                                   ('Gap Extend', 'entry'), ('Output File', 'entry_browse'), (
                                       'Output Format', 'optionmenu',
                                       ["markx0", "markx1", "markx2", "markx3", "markx10", "srspair", "score",
                                        "fasta"]),
                                   ('Tool Website', 'button', 'https://emboss.sourceforge.net/download/'),
                                   ('Version Checked', 'label', 'EMBOSS-6.6.0')],
    },
    "Tree Tools": {
        'RAxML': [('Executable Path', 'entry_browse'), ('Input File', 'entry_browse'),
                  ('Output File', 'entry_browse'), ('Model', 'entry'), ('Seed', 'entry'),
                  ('Tool Website', 'button', 'https://cme.h-its.org/exelixis/web/software/raxml/index.html'),
                  ('Version Checked', 'label', 'version 8.2.12')],
        'FastTree': [('Executable Path', 'entry_browse'), ('Input File', 'entry_browse'),
                     ('Output File', 'entry_browse'), ('Type', 'optionmenu', ["protein", "nucleotide"]),
                     ('Tool Website', 'button', 'http://www.microbesonline.org/fasttree/'),
                     ('Version Checked', 'label', 'Version 2.1.11')],
        'FigTree': [('Executable Path', 'entry_browse'), ('Input File', 'entry_browse'),
                    ('Tool Website', 'button', 'http://tree.bio.ed.ac.uk/software/Figtree/'),
                    ('Version Checked', 'label', 'v1.4.4')
                    ],
    }
}


class Settings(ModalWindow):
    def __init__(self, parent):
        super().__init__(parent,  resizable=True, min_size=(600, 400))
        self.parent = parent
        self.title("Settings")

        self.is_closed = False

        # Create Notebook
        notebook = ttk.Notebook(self)
        notebook.pack(expand=True, fill='both')

        # Create tabs
        self.general_tab = ttk.Frame(notebook)
        self.alignment_tools_tab = ttk.Frame(notebook)
        self.tree_tools_tab = ttk.Frame(notebook)
        self.other_tools_tab = ttk.Frame(notebook)

        notebook.add(self.general_tab, text="General")
        notebook.add(self.alignment_tools_tab, text="Alignment Tools")
        notebook.add(self.tree_tools_tab, text="Tree Tools")
        notebook.add(self.other_tools_tab, text="Other Tools")

        # Load existing settings if any
        self.general_values = load_data(settings_general_lock_path, settings_general_path)
        if self.general_values is None:
            self.general_values = {}

        self.general_entries = {}
        self.create_general_tab()

        # Store input values for each tool (entry, entry_browse and optionmenu values)
        # Dict keys: tool type -> tool -> item
        # Note: "show_tool" stores "widgets" to be able to get most recent values while saving.
        #       "save_current_tool_values" stores "values" as str.
        # Load existing settings if any
        self.tool_values = load_data(settings_tools_lock_path, settings_tools_path)
        if self.tool_values is None:
            self.tool_values = {tool_type: {} for tool_type in tool_params.keys()}

        # Set initial indexes and lists
        self.alignment_tool_index = 0
        self.tree_tool_index = 0
        self.alignment_tool_list = list(tool_params["Alignment Tools"].keys())
        self.tree_tool_list = list(tool_params["Tree Tools"].keys())

        # Display the initial tools in each tab
        self.create_tool_view(self.alignment_tools_tab, "Alignment Tools")
        self.create_tool_view(self.tree_tools_tab, "Tree Tools")

        # Load existing settings if any
        self.othertool_values = load_data(settings_othertools_lock_path, settings_othertools_path)
        if self.othertool_values is None:
            self.othertool_values = []
        self.othertool_frames = []  # for adding frames

        # Create layout for Other Tools tab
        self.create_other_tools_layout()

        # Override the close protocol to save settings
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def on_close(self):
        """Save settings when the Settings window is closed."""
        self.is_closed = True

        # Save General tab
        self.save_general_values()
        save_data(self.general_values, settings_general_lock_path, settings_general_path)

        # Save Alignment Tools and Tree Tools tabs
        self.save_current_tool_values("Alignment Tools")
        self.save_current_tool_values("Tree Tools")
        save_data(self.tool_values, settings_tools_lock_path, settings_tools_path)

        # Save Other Tools tab
        self.save_othertool_values()
        save_data(self.othertool_values, settings_othertools_lock_path, settings_othertools_path)

        # Perform unbind; otherwise, it gives an error.
        self.unbind_mouse_scroll()

        self.destroy()

    def create_general_tab(self):
        validate_cmd = self.general_tab.register(self.validate_int_input)  # Register the validation function with Tkinter

        frame = ttk.Frame(self.general_tab)
        frame.pack(side="top")

        label = ttk.Label(frame, text="Initial Font Size")
        label.grid(row=0, column=0, sticky='w', padx=5, pady=5)
        entry = ttk.Entry(frame,  validate="key", validatecommand=(validate_cmd, "%P"))
        entry.insert(0, self.general_values.get("Font Size", "0"))
        entry.grid(row=0, column=1, sticky='we', padx=5)
        self.general_entries["Font Size"] = entry

        label = ttk.Label(frame, text="Maximum number of operations retained on the undo stack")
        label.grid(row=1, column=0, sticky='w', padx=5, pady=5)
        entry = ttk.Entry(frame,  validate="key", validatecommand=(validate_cmd, "%P"))
        entry.insert(0, self.general_values.get("Undo Limit", "0"))
        entry.grid(row=1, column=1, sticky='we', padx=5)
        self.general_entries["Undo Limit"] = entry

        label = ttk.Label(frame, text="Number of sequences to be cached in index mode")
        label.grid(row=2, column=0, sticky='w', padx=5, pady=5)
        entry = ttk.Entry(frame,  validate="key", validatecommand=(validate_cmd, "%P"))
        entry.insert(0, self.general_values.get("Index Limit", "0"))
        entry.grid(row=2, column=1, sticky='we', padx=5)
        self.general_entries["Index Limit"] = entry

        label = ttk.Label(self.general_tab, text="Note: Options in this tab requires app to be restarted. ", anchor="w")
        label.pack(side="bottom", fill="x")

    # Function to validate integer input
    def validate_int_input(self, value):
        if value == "" or value.isdigit():  # Allows empty input or digits
            return True
        return False

    def save_general_values(self):
        self.general_values = {}

        for name, entry in self.general_entries.items():
            self.general_values[name] = entry.get()

        # check font size
        font_size = check_font_size_range(int(self.general_values["Font Size"]))
        self.general_values["Font Size"] = str(font_size)

    def create_tool_view(self, tab_frame, tool_type): # ("Alignment Tools" and "Tree Tools" tabs)
        # Frame for tool display and navigation buttons
        tool_frame = ttk.Frame(tab_frame)
        tool_frame.pack(fill='both', expand=True)

        # Navigation buttons within the tab
        button_frame = ttk.Frame(tab_frame)
        button_frame.pack(fill='x')
        prev_button = ttk.Button(button_frame, text="Prev", command=lambda: self.show_prev_tool(tool_type))
        prev_button.pack(side='left', padx=5, pady=5)
        next_button = ttk.Button(button_frame, text="Next", command=lambda: self.show_next_tool(tool_type))
        next_button.pack(side='right', padx=5, pady=5)

        # Store button references for enabling/disabling them
        if tool_type == "Alignment Tools":
            self.alignment_prev_button = prev_button
            self.alignment_next_button = next_button
        elif tool_type == "Tree Tools":
            self.tree_prev_button = prev_button
            self.tree_next_button = next_button

        # Store references for each tab to update the tool display
        if tool_type == "Alignment Tools":
            self.alignment_tool_frame = tool_frame
            self.show_tool(tool_type)
        elif tool_type == "Tree Tools":
            self.tree_tool_frame = tool_frame
            self.show_tool(tool_type)

    def show_tool(self, tool_type):
        """
        show tool at current index
        """
        # print("show_tool-------")
        # Clear the existing tool frame for the current tab
        tool_frame, tool_list, tool_index = (None,)*3
        if tool_type == "Alignment Tools":
            tool_frame = self.alignment_tool_frame
            tool_list = self.alignment_tool_list
            tool_index = self.alignment_tool_index
        elif tool_type == "Tree Tools":
            tool_frame = self.tree_tool_frame
            tool_list = self.tree_tool_list
            tool_index = self.tree_tool_index

        for widget in tool_frame.winfo_children():
            widget.destroy()  # only label frame is deleted here

        # Get current tool name and items
        tool_name = tool_list[tool_index]
        tool_items = tool_params[tool_type][tool_name]

        # Display the tool label frame
        label_frame = ttk.LabelFrame(tool_frame, text=tool_name)
        label_frame.pack(fill='x', padx=5, pady=5)

        # Load previous values if available
        tool_values = self.tool_values[tool_type].get(tool_name, {})

        # Create the label frame with the tool's items
        for i, (item_name, item_type, *options) in enumerate(tool_items):
            label = ttk.Label(label_frame, text=item_name)
            label.grid(row=i, column=0, sticky='w', padx=5, pady=5)

            if item_type == 'entry_browse' or item_type == 'entry':
                entry = ttk.Entry(label_frame)
                entry_value = tool_values.get(item_name, "")  # Load saved value or default to empty
                entry.insert(0, entry_value)
                entry.grid(row=i, column=1, sticky='we', padx=5)
                if item_type == 'entry_browse':
                    browse_button = ttk.Button(label_frame, text="Browse",
                                               command=lambda e=entry, item=item_name: self.browse_file(e, item == "Output File"))
                    browse_button.grid(row=i, column=2, padx=5)
                # Store entry reference for saving its value (store widget to get final value while saving)
                tool_values[item_name] = entry

            elif item_type == 'optionmenu':
                var = tk.StringVar()
                var.set(tool_values.get(item_name, options[0][0]))  # Default to the saved or first opt
                # print(item_name, var.get(), options)
                option_menu = tk.OptionMenu(label_frame, var, *options[0])
                option_menu.grid(row=i, column=1, sticky='w', padx=5)
                # Store variable reference for saving its value (store widget to get final value while saving)
                tool_values[item_name] = var

            elif item_type == 'button':
                browse_button = ttk.Button(label_frame, text="Open",
                                           command=lambda url=options[0]: webbrowser.open(url))   # option's curr value
                browse_button.grid(row=i, column=1, sticky='w', padx=5)

            elif item_type == 'label':
                label = ttk.Label(label_frame, text=options[0])
                label.grid(row=i, column=1, sticky='w', padx=5, pady=5)

            label_frame.columnconfigure(1, weight=1)  # Make entry columns expandable

        # Save updated tool values reference
        self.tool_values[tool_type][tool_name] = tool_values

        # Update button states
        self.update_navigation_buttons(tool_type)

    def update_navigation_buttons(self, tool_type):
        # Disable "Prev" on the first tool, "Next" on the last tool in each tab
        if tool_type == "Alignment Tools":
            self.alignment_prev_button.config(state="normal" if self.alignment_tool_index > 0 else "disabled")
            self.alignment_next_button.config(state="normal" if self.alignment_tool_index < len(self.alignment_tool_list) - 1 else "disabled")
        elif tool_type == "Tree Tools":
            self.tree_prev_button.config(state="normal" if self.tree_tool_index > 0 else "disabled")
            self.tree_next_button.config(state="normal" if self.tree_tool_index < len(self.tree_tool_list) - 1 else "disabled")

    def show_next_tool(self, tool_type):
        self.save_current_tool_values(tool_type)
        if tool_type == "Alignment Tools" and self.alignment_tool_index < len(self.alignment_tool_list) - 1:
            self.alignment_tool_index += 1
            self.show_tool(tool_type)
        elif tool_type == "Tree Tools" and self.tree_tool_index < len(self.tree_tool_list) - 1:
            self.tree_tool_index += 1
            self.show_tool(tool_type)

    def show_prev_tool(self, tool_type):
        self.save_current_tool_values(tool_type)
        if tool_type == "Alignment Tools" and self.alignment_tool_index > 0:
            self.alignment_tool_index -= 1
            self.show_tool(tool_type)
        elif tool_type == "Tree Tools" and self.tree_tool_index > 0:
            self.tree_tool_index -= 1
            self.show_tool(tool_type)

    def save_current_tool_values(self, tool_type):
        # print("save_current_tool_values-------")
        # Save the values in entries and optionmenus for the currently displayed tool (called when user clicks prev or next)
        tool_list, tool_index = None, None
        if tool_type == "Alignment Tools":
            tool_list = self.alignment_tool_list
            tool_index = self.alignment_tool_index
        elif tool_type == "Tree Tools":
            tool_list = self.tree_tool_list
            tool_index = self.tree_tool_index

        tool_name = tool_list[tool_index]
        tool_values = self.tool_values[tool_type][tool_name]

        # for item_name, widget in tool_values.items():
        #     # print(item_name, widget.get())
        #     if isinstance(widget, tk.Entry):
        #         # Save entry text
        #         self.tool_values[tool_type][tool_name][item_name] = widget.get()  # What if we use tool_values directly?
        #     elif isinstance(widget, tk.StringVar):
        #         # Save option menu selection
        #         self.tool_values[tool_type][tool_name][item_name] = widget.get()

        for item_name, widget in list(tool_values.items()):  # Create a list copy of items to change while iterating
            if isinstance(widget, tk.Entry):
                tool_values[item_name] = widget.get()
            elif isinstance(widget, tk.StringVar):
                tool_values[item_name] = widget.get()

    def browse_file(self, entry_widget, is_save):
        if is_save:
            filepath = filedialog.asksaveasfilename()
        else:
            filepath = filedialog.askopenfilename()
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, filepath)

    def create_other_tools_layout(self):
        # Scrollable frame setup
        self.canvas = tk.Canvas(self.other_tools_tab, highlightthickness=0, bd=0)
        self.canvas.pack(side="left", fill="both", expand=True)

        self.scrollable_frame = ttk.Frame(self.other_tools_tab)
        self.last_frame_id = None
        self.last_width = None

        # Initial tool row - cannot be deleted
        if len(self.othertool_values) > 0:
            self.add_first_tool_row(*self.othertool_values[0])
        else:
            self.add_first_tool_row()

        # other tool rows
        for tool_name, tool_command in self.othertool_values[1:]:
            self.add_tool_row(tool_name, tool_command)

        # Add New button
        add_button = ttk.Button(self.scrollable_frame, text="Add New", command=self.add_new_tool_row)
        add_button.pack(side="bottom", pady=4, padx=10)

        # Scrollbar
        self.v_scrollbar = ttk.Scrollbar(self.other_tools_tab, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.v_scrollbar.set)

        self.bind_mouse_scroll()

        # Pack the canvas and scrollbar
        self.v_scrollbar.pack(side="right", fill="y")

        self.other_tools_tab.bind('<Configure>', self.update_size)

    def update_size(self, event):
        # height = max(self.other_tools_height, self.other_tools_tab.winfo_height())
        # height = self.other_tools_height
        # print(self.scrollable_frame.winfo_height(), self.other_tools_tab.winfo_height())
        if self.is_closed:
            return

        # enable/disable scroll
        list_height = self.scrollable_frame.winfo_height()
        window_height = self.other_tools_tab.winfo_height()
        if list_height == 1 or list_height > window_height:  # list_height is 1 at startup
            self.bind_mouse_scroll()
            self.v_scrollbar.pack(side="right", fill="y")
        else:
            self.unbind_mouse_scroll()
            self.v_scrollbar.pack_forget()

        # recreate scrollable_frame if width changed
        width = self.other_tools_tab.winfo_width()
        if self.last_width != width:
            # delete old one
            if self.last_frame_id is not None:
                self.canvas.delete(self.last_frame_id)

            # add new one
            self.last_frame_id = self.canvas.create_window(
                                 (0, 0),
                                 window=self.scrollable_frame,
                                 anchor='nw',
                                 width=width,
            )
        self.last_width = width

        try:
            self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        except: # gives exception while window is being closed
            pass

    def bind_mouse_scroll(self):
        self.canvas.bind_all('<MouseWheel>',
                             lambda event: self.canvas.yview_scroll(-int(event.delta / 1), "units"))

    def unbind_mouse_scroll(self):
        self.canvas.unbind_all('<MouseWheel>')

    def add_first_tool_row(self, tool_name="", tool_command=""):
        """Adds a new row for a tool."""
        frame = ttk.Frame(self.scrollable_frame)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=0)
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(2, weight=0)

        # Header row
        ttk.Label(frame, text="Tool Name").grid(row=0, column=0, sticky='w',padx=5, pady=5)
        ttk.Label(frame, text="Tool Command").grid(row=0, column=1, sticky='w',padx=5, pady=5)

        # Tool name entry
        tool_name_entry = ttk.Entry(frame, width=10)
        tool_name_entry.grid(row=1, column=0,  padx=5, pady=5)
        tool_name_entry.insert(0, tool_name)

        # Tool command entry
        tool_command_entry = ttk.Entry(frame)
        tool_command_entry.grid(row=1, column=1, sticky='we', padx=5, pady=5)
        tool_command_entry.insert(0, tool_command)

        # Remove button
        remove_button = ttk.Button(frame, text="-")
        remove_button.grid(row=1, column=2, padx=10, pady=5)

        frame.pack(fill='both', pady=4, padx=10)

        # Append frame and entries to the list
        self.othertool_frames.append((frame, tool_name_entry, tool_command_entry))

    def add_tool_row(self, tool_name="", tool_command=""):
        """Adds a new row for a tool."""
        frame = ttk.Frame(self.scrollable_frame)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=0)
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(2, weight=0)

        # Tool name entry
        tool_name_entry = ttk.Entry(frame, width=10)
        tool_name_entry.grid(row=0, column=0,  padx=5, pady=5)
        tool_name_entry.insert(0, tool_name)

        # Tool command entry
        tool_command_entry = ttk.Entry(frame)
        tool_command_entry.grid(row=0, column=1, sticky='we', padx=5, pady=5)
        tool_command_entry.insert(0, tool_command)

        # Remove button
        remove_button = ttk.Button(frame, text="Remove", command=lambda: self.remove_tool_row(frame))
        remove_button.grid(row=0, column=2, padx=10, pady=5)

        frame.pack(fill='both', pady=4, padx=10)

        # Append frame and entries to the list
        self.othertool_frames.append((frame, tool_name_entry, tool_command_entry))

    def add_new_tool_row(self):
        self.add_tool_row()
        self.scrollable_frame.update_idletasks()
        self.update_size(None)

    def remove_tool_row(self, frame):
        # Find and remove frame from list
        self.othertool_frames = [(f, e1, e2) for f, e1, e2 in self.othertool_frames if f != frame]
        frame.destroy()
        self.scrollable_frame.update_idletasks()
        self.update_size(None)

    def save_othertool_values(self):
        self.othertool_values = []

        for frame, entry1, entry2 in self.othertool_frames:
            self.othertool_values.append((entry1.get(), entry2.get()))


def load_data(lock_path, file_path):
    """
    Load tool values from settings file if it exists.
    We used FileLock with 0.1 sec blocking to prevent multiple program instances from conflicting
    """
    if os.path.exists(file_path):
        try:
            with FileLock(lock_path, timeout=0.1):
                with open(file_path, "r") as file:
                    data = json.load(file)
                    # print("Data loaded successfully:", file_path)
                    return data
        except Timeout:
            print("File is currently locked by another instance. Could not load.")
    return None


def save_data(data, lock_path, file_path):
    """
    Save self.tool_values to settings file as JSON.
    We used FileLock with 0.1 sec blocking to prevent multiple program instances from conflicting
    """
    try:
        with FileLock(lock_path, timeout=0.1):  # Timeout after 0.1 seconds if unable to acquire lock
            with open(file_path, "w") as file:
                json.dump(data, file)
                # print("Tool values saved successfully.")
    except Timeout:
        print("File is currently locked by another instance. Could not save tool values.")
