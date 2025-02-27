import shlex
import subprocess
import threading
import webbrowser
from time import sleep
from tkinter import messagebox, ttk, filedialog

from filelock import FileLock, Timeout

from utils.model_window import ModalWindow
from utils.windows import TKMessageBox
from utils.settings import tool_params, settings_tools_lock_path, save_data, settings_othertools_lock_path, \
    settings_othertools_path, load_data
import os
import json
import tkinter as tk
from tkinter import ttk, filedialog
from appdirs import AppDirs
from utils.settings import settings_tools_path

OUTPUT_TO_FILE = 1  # tools that write aligned seqs to stdout.
OUTPUT_TO_STDOUT = 2  # tools that write aligned seqs to stdout.

prank_format_dict = {
    "IG/Stanford": 1,
    "Pearson/Fasta": 8,
    "GenBank/GB": 2,
    "Phylip3.2": 11,
    "NBRF": 3,
    "Phylip": 12,
    "EMBL": 4,
    "PIR/CODATA": 14,
    "DNAStrider": 6,
    "MSF": 15,
    "Fitch": 7,
    "PAUP/NEXUS": 17
}


class ToolProcessManager(ModalWindow):
    # Buffered I/O: When a process writes to stdout or stderr, the data is not necessarily sent to its destination
    # (like a terminal, file, or pipe) immediately. Instead, the data is often temporarily stored in a buffer
    # (a memory area) before being written to the final output.
    # By default, stdout is line-buffered when writing to a terminal. This means the output is flushed
    # (sent to the destination) when a newline character (\n) is encountered. If no newline is written,
    # the data stays in the buffer until it is either flushed manually or the buffer fills up.
    # In contrast, stderr is usually unbuffered, meaning data written to stderr is immediately flushed and becomes
    # available for reading right away.
    # "Flushing" refers to forcing the data in the buffer to be written to the destination (e.g., a file or pipe).
    # A process can explicitly flush its buffers by calling flush(), or the buffer may be automatically flushed when:
    #
    # - A newline character is written (in line-buffered mode).
    # - The buffer becomes full.
    # - The process terminates.
    def __init__(self, parent, tool_type: int, command:list, output_file_path=None):
        self.parent = parent  # main window
        super().__init__(self.parent, min_size=(200, 150))
        self.tool_type = tool_type
        self.command = command
        self.output_file_path = output_file_path

        self.title("Progress")

        # Status frame
        status_frame = tk.Frame(self)
        status_frame.pack(side=tk.BOTTOM, fill=tk.BOTH)
        self.status_var = tk.StringVar(value="Running...")
        status_label = tk.Label(status_frame, font=("Arial", 12), textvariable=self.status_var)
        status_label.pack()

        # Progress frame
        progress_frame = tk.Frame(self)
        progress_frame.pack(side=tk.TOP, expand=True, fill=tk.BOTH)

        tw_frame = tk.Frame(progress_frame)
        tw_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # tw
        self.text_widget = tk.Text(tw_frame, wrap=tk.NONE, state="disabled",
                                   font=("Courier", 16), bd=0, highlightthickness=0)
        self.text_widget.pack(expand=True, fill=tk.BOTH)

        scrollbar_x = tk.Scrollbar(tw_frame, orient=tk.HORIZONTAL, command=self.text_widget.xview)
        scrollbar_y = tk.Scrollbar(progress_frame, orient=tk.VERTICAL, command=self.text_widget.yview)
        self.text_widget.config(yscrollcommand=scrollbar_y.set, xscrollcommand=scrollbar_x.set)

        scrollbar_x.place(relx=0, rely=1, relwidth=1, anchor='sw')
        scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)

        self.tool_process = None
        self.running = False
        self.output_file = None  # file for tools
        self.window_closed = False

        self.protocol("WM_DELETE_WINDOW", self.on_close)

        # start tool after creating progress window
        self.append_text("Command:\n" + " ".join(self.command) + "\n\n")
        self.start_process()

    def start_process(self):
        self.running = True
        # Open a file to redirect stderr
        try:
            if self.tool_type == OUTPUT_TO_FILE:
                self.tool_process = subprocess.Popen(self.command,
                                                     stdout=subprocess.PIPE,  # progress --> text widget
                                                     stderr=subprocess.PIPE,  # error --> text widget
                                                     text=True)
            elif self.tool_type == OUTPUT_TO_STDOUT:
                self.output_file = open(self.output_file_path, 'w')
                self.tool_process = subprocess.Popen(self.command,
                                                     stdout=self.output_file,  # output --> file
                                                     stderr=subprocess.PIPE,  # progress and/or error --> text widget
                                                     text=True)
            else:
                raise Exception("Not valid tool type.")
        except Exception as e:
            self.running = False
            self.append_text(str(e) + "\n")
            self.status_var.set(f"Stopped.")
            return

        # Start a thread to read the progress from stdout
        if self.tool_type == OUTPUT_TO_FILE:
            threading.Thread(target=self.read_std_out, daemon=True).start()

        # Start a thread to read the progress(OUTPUT_TO_STDOUT) or/and error(OUTPUT_TO_FILE) from stderr
        threading.Thread(target=self.read_std_err, daemon=True).start()

    def read_std_out(self):
        """
        Run in OUTPUT_TO_FILE mode
        Read stdout PIPE, transfer into text widget.
        :return:
        """
        pipe = self.tool_process.stdout

        # read pipe until process ends
        while self.running:
            line = pipe.readline()  # blocked if no available data
            if line:
                self.master.after(0, self.append_text, line)
                # print("line:", line)
            elif self.tool_process.poll() is not None:  # Process has finished
                self.running = False
                break

        # Check for any remaining output after the process has ended. (possible due to buffered IO)
        # print("remaining")
        for line in pipe:
            self.master.after(0, self.append_text, line)
            # print("line:", line)

        # Clean up: close the stdout
        if pipe:
            pipe.close()

        print("read_std_out finished")

    def read_std_err(self):
        """
         Run in BOTH mode
         Read stdout PIPE, transfer into text widget.
         :return:
         """
        pipe = self.tool_process.stderr

        while self.running:
            line = pipe.readline()  # blocked if no available data
            if line:
                self.master.after(0, self.append_text, line)
            elif self.tool_process.poll() is not None:  # Process has finished
                self.running = False
                break

        # Check for any remaining output after the process has ended. (due to buffered IO)
        for line in pipe:
            self.master.after(0, self.append_text, line)

        # Clean up: close the stderr
        if pipe:
            pipe.close()

        # Clean up: close the file
        if self.tool_type == OUTPUT_TO_STDOUT and self.output_file:
            self.output_file.close()

        # Print "Completed" if tool has finished successfully
        if not self.window_closed:
            if self.tool_process.returncode == 0:
                self.status_var.set("Completed successfully.")
            else:
                self.status_var.set(f"Stopped with error code {self.tool_process.returncode}.")

        print("read_std_err finished")

    def append_text(self, text):
        if not self.window_closed:
            self.text_widget.configure(state="normal")
            self.text_widget.insert(tk.END, text)
            self.text_widget.configure(state="disabled")
            self.text_widget.see(tk.END)   # Scroll to the end

    def on_close(self):
        if self.running:
            quit_tool = TKMessageBox("Warning", "The tool is currently running. Do you want to quit it now ?", messagebox.askyesno, self).show()
            if quit_tool:
                if self.tool_process:
                    self.tool_process.terminate()
            else:
                return  # kapatmaktan vazgeçti
        self.window_closed = True
        self.destroy()

    def disable_menu(self):
        pass

    def enable_menu(self):
        pass


class ToolWindow(ModalWindow):
    def __init__(self, parent, tool_type, tool_name, os_type=None):
        super().__init__(parent, resizable=False)
        self.title(f"{tool_name}")
        self.tool_type = tool_type
        self.tool_name = tool_name
        self.os_type = os_type

        # Label frame for the specific tool
        label_frame = ttk.LabelFrame(self, text=tool_name)
        label_frame.pack(fill='both', expand=True, padx=10, pady=10)

        # Get tool items for the specific tool
        tool_items = tool_params[tool_type][tool_name]

        # Initialize tool values
        self.all_values = load_data(settings_tools_lock_path, settings_tools_path)
        if self.all_values is None:
            self.all_values = {}

        self.tool_values = self.all_values.get(tool_type, {}).get(tool_name, {})

        # Populate the label frame with tool settings
        for i, (item_name, item_type, *options) in enumerate(tool_items):
            label = ttk.Label(label_frame, text=item_name)
            label.grid(row=i, column=0, sticky='w', padx=5, pady=5)

            if item_type == 'entry_browse' or item_type == 'entry':
                entry = ttk.Entry(label_frame)
                entry_value = self.tool_values.get(item_name, "")  # Load saved value or default to empty
                entry.insert(0, entry_value)
                entry.grid(row=i, column=1, sticky='we', padx=5)
                if item_type == 'entry_browse':
                    # print(item_name)
                    browse_button = ttk.Button(label_frame, text="Browse",
                                               command=lambda e=entry, item=item_name: self.browse_file(e, item == "Output File"))
                    browse_button.grid(row=i, column=2, padx=5)
                self.tool_values[item_name] = entry

            elif item_type == 'optionmenu':
                var = tk.StringVar()
                var.set(self.tool_values.get(item_name, options[0][0]))  # Default to the first option
                option_menu = tk.OptionMenu(label_frame, var, *options[0])
                option_menu.grid(row=i, column=1, sticky='w', padx=5)
                self.tool_values[item_name] = var

            elif item_type == 'button':
                browse_button = ttk.Button(label_frame, text="Open",
                                           command=lambda url=options[0]: webbrowser.open(url))  # option's curr value
                browse_button.grid(row=i, column=1, sticky='w', padx=5)

            elif item_type == 'label':
                label = ttk.Label(label_frame, text=options[0])
                label.grid(row=i, column=1, sticky='w', padx=5, pady=5)

            label_frame.columnconfigure(1, weight=1)  # Make entry columns expandable

        # Run button centered at the bottom
        run_button = ttk.Button(self, text="Run", command=self.run_tool)
        run_button.pack(pady=(10, 20))  # Add padding to center the button at the bottom

        # Override the close protocol to save settings
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def on_close(self):
        # Save tool values
        self.save_tool_values()
        self.destroy()

    def save_tool_values(self):
        """
        Save the current tool's values to the main settings file.
        We used FileLock with 0.1 sec blocking to prevent multiple program instances from conflicting
        """

        # Update the specific tool's values
        if self.tool_type not in self.all_values:
            self.all_values[self.tool_type] = {}
        self.all_values[self.tool_type][self.tool_name] = {
            name: widget.get() for name, widget in self.tool_values.items()
        }

        save_data(self.all_values, settings_tools_lock_path, settings_tools_path)

    def run_tool(self):
        params = [widget.get() for widget in self.tool_values.values()]

        if self.tool_name == "Muscle":
            cmd = muscle_command_line(*params)
            ToolProcessManager(self.parent, OUTPUT_TO_FILE, cmd)
        elif self.tool_name == "ClustalOmega":
            cmd = clustalomega_command_line(*params)
            ToolProcessManager(self.parent, OUTPUT_TO_FILE, cmd)
        elif self.tool_name == "Prank":
            params[-1] = prank_format_dict[params[-1]]  # convert format to number
            cmd = prank_command_line(*params)
            ToolProcessManager(self.parent, OUTPUT_TO_FILE, cmd)
        elif self.tool_name == "Mafft":
            if params[-1] == "phylip":
                phylipout = True
                clustalout = False
            elif params[-1] == "clustal":
                phylipout = False
                clustalout = True
            else:
                phylipout = False
                clustalout = False
            cmd = mafft_command_line(*params[:-2], phylipout, clustalout)   # exclude outfile and outformat
            ToolProcessManager(self.parent, OUTPUT_TO_STDOUT, cmd, params[-2])
        elif self.tool_name == "Probcons":
            if params[-1] == "ClustalW":
                clustalw = True
            else:
                clustalw = False
            cmd = probcons_command_line(*params[:-2], clustalw)
            ToolProcessManager(self.parent, OUTPUT_TO_STDOUT, cmd, params[-2])
        elif self.tool_name == "MSAProbs":
            if params[-1] == "CLUSTALW":
                clustalw = True
            else:
                clustalw = False
            cmd = msaprobs_command_line(*params[:-1], clustalw)
            ToolProcessManager(self.parent, OUTPUT_TO_FILE, cmd)
        elif self.tool_name == "Needleman-Wunsch (Global)":
            cmd = needle_command_line(*params)
            ToolProcessManager(self.parent, OUTPUT_TO_FILE, cmd)
        elif self.tool_name == "Smith-Waterman (Local)":
            cmd = water_command_line(*params)
            ToolProcessManager(self.parent, OUTPUT_TO_FILE, cmd)
        elif self.tool_name == "RAxML":
            out_file_dir = os.path.dirname(params[2])
            out_file_name = os.path.basename(params[2])
            cmd = raxml_command_line(params[0], params[1], out_file_dir, out_file_name, params[3], params[4])
            ToolProcessManager(self.parent, OUTPUT_TO_FILE, cmd)
        elif self.tool_name == "FastTree":
            nucleotide = True if params[3] == "nucleotide" else False
            cmd = fasttree_command_line(params[0], params[1], nucleotide)
            ToolProcessManager(self.parent, OUTPUT_TO_STDOUT, cmd, params[2])
        elif self.tool_name == "FigTree":
            if self.os_type is None:
                raise Exception("OS type not defined !")
            cmd = figtree_command_line(*params, self.os_type)
            subprocess.Popen(cmd)  # exe, main script, open file mode

        self.on_close()

    def browse_file(self, entry_widget, is_save=False):
        if is_save:
            filepath = filedialog.asksaveasfilename()
        else:
            filepath = filedialog.askopenfilename()
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, filepath)


class OtherToolWindow(ModalWindow):
    def __init__(self, parent, tool_index, tool_name, tool_command, other_tools_data):
        super().__init__(parent, resizable=False)
        self.title(f"{tool_name}")

        self.tool_index = tool_index
        self.other_tools_data = other_tools_data

        # Main frame setup
        frame = ttk.Frame(self)
        frame.pack(fill='both', pady=4, padx=10, expand=True)

        # Configure the frame grid
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=0)
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(2, weight=0)

        # Header row
        ttk.Label(frame, text="Tool Name").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        ttk.Label(frame, text="Tool Command").grid(row=0, column=1, sticky='w', padx=5, pady=5)

        # Tool name entry
        self.tool_name_entry = ttk.Entry(frame, width=10)
        self.tool_name_entry.grid(row=1, column=0, padx=5, pady=5)
        self.tool_name_entry.insert(0, tool_name)

        # Tool command entry
        self.tool_command_entry = ttk.Entry(frame, width=30)
        self.tool_command_entry.grid(row=1, column=1, sticky='we', padx=5, pady=5)
        self.tool_command_entry.insert(0, tool_command)

        # Run button
        run_button = ttk.Button(self, text="Run", command=self.run_tool)
        run_button.pack(side='bottom', pady=10)

        # Override the close protocol to save settings
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def on_close(self):
        # Save tool values
        self.save_tool_values()
        self.destroy()

    def run_tool(self):
        tool_name = self.tool_name_entry.get()
        tool_command = self.tool_command_entry.get()
        print(f"Running tool: {tool_name} with command: {tool_command}")

        # run
        command_list = shlex.split(tool_command)
        ToolProcessManager(self.parent, OUTPUT_TO_FILE, command_list)

        self.on_close()

    def save_tool_values(self):
        tool_name = self.tool_name_entry.get()
        tool_command = self.tool_command_entry.get()
        self.other_tools_data[self.tool_index] = (tool_name, tool_command)  # update tool
        save_data(self.other_tools_data, settings_othertools_lock_path, settings_othertools_path)


def muscle_command_line(exe, in_file, out_file, out_format):
    cmd = [exe, "-align", in_file, "-output", out_file]
    return cmd


def clustalomega_command_line(exe, in_file, out_file, out_format):
    # --force: Force file overwriting.
    cmd = [exe, "-i", in_file, "-o", out_file, "--outfmt", out_format, "-v", "--force"]
    return cmd


def prank_command_line(exe, in_file, out_file, out_format):
    cmd = [exe, f"-d={in_file}", f"-o={out_file}", f"-f={out_format}"]
    return cmd


def mafft_command_line(exe, in_file, phylipout, clustalout):
    if phylipout:
        cmd = [exe, "--phylipout", in_file]
    elif clustalout:
        cmd = [exe, "--clustalout", in_file]
    else:
        cmd = [exe, in_file]
    return cmd


def probcons_command_line(exe, in_file, clustalw):
    if clustalw:
        cmd = [exe, "-clustalw", "--verbose", in_file]
    else:
        cmd = [exe, "--verbose", in_file]
    return cmd


def msaprobs_command_line(exe, in_file, out_file, clustalw):
    if clustalw:
        cmd = [exe, in_file, "-o", out_file, "-v", "-clustalw"]
    else:
        cmd = [exe, in_file, "-o", out_file, "-v"]
    return cmd


def needle_command_line(exe, in_file_1, in_file_2, gapopen, gapextend, out_file, out_format):
    cmd = [exe, "-verbose", f"-outfile={out_file}", f"-asequence={in_file_1}", f"-bsequence={in_file_2}", f"-gapopen={gapopen}", f"-gapextend={gapextend}", f"-aformat={out_format}"]
    return cmd


def water_command_line(exe, in_file_1, in_file_2, gapopen, gapextend, out_file, out_format):
    cmd = [exe, "-verbose", f"-outfile={out_file}", f"-asequence={in_file_1}", f"-bsequence={in_file_2}", f"-gapopen={gapopen}", f"-gapextend={gapextend}", f"-aformat={out_format}"]
    return cmd


def raxml_command_line(exe, in_file, out_file_dir, out_file_name, model, seed):
    cmd = [exe, "-s", in_file, "-w", out_file_dir, "-n", out_file_name, "-m", model, "-p", seed]
    return cmd


def fasttree_command_line(exe, in_file, nucleotide):
    if nucleotide:
        cmd = [exe, "-nt", in_file]
    else:
        cmd = [exe, in_file]
    return cmd


def figtree_command_line(exe, in_file, os_type):
    if os_type == "mac":
        cmd = ["open", "-a", exe, in_file]  # not works
    elif os_type == "windows":
        cmd = [exe, in_file]  # works
    else:
        cmd = [exe, in_file]
    return cmd



