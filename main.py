import platform
import tkinter as tk
from tkinter import messagebox, ttk, filedialog
import subprocess
import sys

from common.configs import MIN_FONT_SIZE, MAX_FONT_SIZE, DEFAULT_FONT_SIZE
from common.data import Data
from common.shortcuts import shortcuts
from external_tools.tools import muscle_command_line, ToolProcessManager, OUTPUT_TO_FILE, \
    clustalomega_command_line, prank_command_line, needle_command_line, water_command_line, OUTPUT_TO_STDOUT, \
    mafft_command_line, probcons_command_line, msaprobs_command_line, fasttree_command_line, figtree_command_line, \
    raxml_command_line, ToolWindow, OtherToolWindow
from tree.tree_convert import ConvertTree
from Bio import motifs as bio_motif
from Bio import SeqIO, AlignIO, Phylo
from sequence_alignment.alignment_file_dialog import AligmentOpenFileDialog, AlignmentSaveFileDialog
from sequence_alignment.alignment_format_converter import ConvertAlignmentFormat
from sequence_editor.sequence_file_converter import ConvertSequenceFormat
from utils.model_window import ModalWindow
from utils.settings import Settings, load_data, settings_othertools_lock_path, settings_othertools_path, \
    settings_general_lock_path, settings_general_path
from utils.windows import TKMessageBox, MessageBox
from viewer import sequence_viewer
from viewer.sequence_viewer import SequenceViewer
from sequence_editor.sequence_reader import SequenceReader, SequenceIndexReader
from sequence_editor.sequence_file_dialog import SequenceOpenFileDialog, SequenceOpenIndexFileDialog, SequenceSaveFileDialog
from sequence_alignment.alignment_reader import AlignmentReader
from utils.utils import toggle_menu_on_off, toggle_menu_on_off_all
import os
from sequence_alignment.alignment_selector import SelectAlignment
from motif.motif_finder import MotifFinder, SelectEditMode
from motif.motif import MotifSearch
from motif.motif_reader import MotifFileDialog
from sequence_editor.sequence_writer import SequenceWriter, SequenceIndexWriter
from sequence_alignment.alignment_writer import AlignmentWriter
import threading
from memory_profiler import profile


OPEN_SEQ_FILE = "OPEN_SEQ_FILE"
OPEN_SEQ_FILE_INDEX = "OPEN_SEQ_FILE_INDEX"
OPEN_ALIGN_FILE = "OPEN_ALIGN_FILE"

DEFAULT_INDEX_LIMIT = 0


class App(tk.Tk):
    # @profile
    def __init__(self, title, size):

        # main setup
        super().__init__()
        self.title(title)
        self.window_width = size[0]
        self.window_height = size[1]

        # Initialize
        self.file_path = None
        self.data = None
        self.unsaved_file = False
        self.read_write_lock = threading.Lock()
        self.message_box = None
        self.save_error = None
        self.motifs = None
        self.motif_finder_open = False
        self.os_type = "mac" if platform.system() == "Darwin" else "windows" if platform.system() == "Windows" else "linux"
        self.load_settings()

        # create viewer
        self.viewer = SequenceViewer(self, self.general_settings)

        # menu
        self.menu_bar = None
        self.menu_names = ["File", "Edit", "Select", "View", "Sequence", "Alignment", "Annotation", "Tree", "Motif", "Other Tools"]
        self.file_menu, self.edit_menu, self.sequence_menu, self.alignment_menu, self.tree_menu, self.motif_menu, self.other_tools_menu = (None,) * 7
        self.create_menu()  # all ENABLED
        self.disable_menu_default()  # DISABLE for default
        self.menu_enabled = True  # False if all DISABLED (when a child window appears)

        # other tools menu
        # self.update_other_tools_menu()
        self.other_tools_data = None

        # Bind keyboard shortcuts
        self.bind_keyboard_shortcuts()

        # Set the protocol to intercept the close event
        self.protocol("WM_DELETE_WINDOW", self.on_closing)

        # check if bundled
        self.is_bundled = self.check_if_bundled()

        # check open mode (if a this is a new instance created with "Open ..." button )
        self.check_open_mode()

    def check_if_bundled(self):
        # print(sys.executable, sys.argv)

        if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
            print('running in a PyInstaller bundle')
            return True
        else:
            print('running in a normal Python process')
            return False

    def load_settings(self):
        general_settings = load_data(settings_general_lock_path, settings_general_path)
        self.general_settings = {} if general_settings is None else general_settings

    def update_other_tools_menu(self):
        if self.menu_enabled:
            self.other_tools_menu.delete(0, tk.END)  # delete all

            # update menu after loading settings
            self.other_tools_data = load_data(settings_othertools_lock_path, settings_othertools_path)
            if self.other_tools_data is not None:
                for i, (tool_name, tool_command) in enumerate(self.other_tools_data):
                    if tool_name != "":
                        self.other_tools_menu.add_command(label=tool_name,
                                                          command=lambda index=i, name=tool_name, command=tool_command:
                                                          self.open_other_tool(index, name, command))

    def open_other_tool(self, tool_index, tool_name, tool_command):
        print(tool_index, tool_name, tool_command)
        OtherToolWindow(self, tool_index, tool_name, tool_command, self.other_tools_data)

    def on_closing(self):
        # save or discard unsaved file
        print("closing app...")
        choice = self.check_unsaved_file()
        if choice is None:
            return  # cancelled
        elif choice:
            self.save_file()  # save.
        else:
            self.destroy()  # no save and close.

    def check_open_mode(self):
        x, y = 100, 100  # default win offsets

        if len(sys.argv) > 1:  # opened by another instance
            window_pos = sys.argv[1].split("=")[1]
            x, y = map(int, window_pos.split(","))

        self.geometry(f'{self.window_width}x{self.window_height}+{x}+{y}')
        self.minsize(self.window_width, self.window_height)

        open_mode = None
        if len(sys.argv) > 2:  # open mode has been set
            open_mode = sys.argv[2]

        if open_mode is not None:
            if open_mode == OPEN_SEQ_FILE:
                self.open_seq_file()
            elif open_mode == OPEN_SEQ_FILE_INDEX:
                self.open_seq_file_index()
            elif open_mode == OPEN_ALIGN_FILE:
                self.open_align_file()

    def create_menu(self):
        # Create menu bar
        self.menu_bar = tk.Menu()
        self.config(menu=self.menu_bar)

        # Create "File" menu
        self.file_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.file_menu.add_command(label="New", command= self.open_new_instance, accelerator=shortcuts["New"][self.os_type])
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Open Sequence File", command=self.open_seq_file, accelerator=shortcuts["Open Sequence File"][self.os_type])
        self.file_menu.add_command(label="Open Sequence File as Index", command=self.open_seq_file_index, accelerator=shortcuts["Open Sequence File as Index"][self.os_type])
        self.file_menu.add_command(label="Open Alignment File", command=self.open_align_file, accelerator=shortcuts["Open Alignment File"][self.os_type])
        self.file_menu.add_command(label="Save File", command=self.save_file, accelerator=shortcuts["Save File"][self.os_type])
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Convert Sequence File Format", command=self.convert_seq_format)
        self.file_menu.add_command(label="Convert Alignment File Format", command=self.convert_align_format)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Settings", command=self.open_settings)
        self.menu_bar.add_cascade(label=self.menu_names[0], menu=self.file_menu)

        # Create a "Edit" menu
        self.edit_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.edit_menu.add_checkbutton(label="Editor Mode (Keyboard)", variable=self.viewer.editor_mode_var,
                                       command=self.viewer.toggle_editor_mode,
                                       accelerator=shortcuts['Editor Mode'][self.os_type])
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Undo", command=self.viewer.undo, accelerator=shortcuts["Undo"][self.os_type])
        self.edit_menu.add_command(label="Redo", command=self.viewer.redo, accelerator=shortcuts["Redo"][self.os_type])
        self.edit_menu.add_command(label="Save Current State", command=self.viewer.save_current_state)

        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Copy Selection as FASTA", command=self.viewer.copy_sel_as_fasta, accelerator=shortcuts["Copy Selection as FASTA"][self.os_type])
        self.edit_menu.add_command(label="Copy Selection as Text", command=self.viewer.copy_sel_as_text, accelerator=shortcuts["Copy Selection as Text"][self.os_type])
        self.edit_menu.add_command(label="Paste FASTA", command=self.viewer.paste_fasta, accelerator=shortcuts["Paste FASTA"][self.os_type])
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Delete Selected", command=self.viewer.menu_delete_selected, accelerator=shortcuts["Delete Selected"][self.os_type])
        self.edit_menu.add_command(label="Delete Selected Letters", command=self.viewer.menu_delete_selected_letters, accelerator=shortcuts["Delete Selected Letters"][self.os_type])
        self.edit_menu.add_command(label="Insert Unknown Letter", command=self.viewer.insert_unknown_letter, accelerator=shortcuts["Insert Unknown Letter"][self.os_type])
        self.menu_bar.add_cascade(label=self.menu_names[1], menu=self.edit_menu)

        # Create a "Select" menu
        self.select_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.select_menu.add_command(
            label="De-select",
            command=self.viewer.de_select,
            accelerator=shortcuts["De-select"][self.os_type]
        )
        self.select_menu.add_command(
            label="Select All",
            command=self.viewer.select_all,
            accelerator=shortcuts["Select All"][self.os_type]
        )
        self.select_menu.add_command(
            label="Expand Selection Right",
            command=self.viewer.expand_sel_right,
            accelerator=shortcuts["Expand Selection Right"][self.os_type]
        )
        self.select_menu.add_command(
            label="Expand Selection Left",
            command=self.viewer.expand_sel_left,
            accelerator=shortcuts["Expand Selection Left"][self.os_type]
        )
        self.select_menu.add_command(
            label="Expand Selection Up",
            command=self.viewer.expand_sel_up,
            accelerator=shortcuts["Expand Selection Up"][self.os_type]
        )
        self.select_menu.add_command(
            label="Expand Selection Down",
            command=self.viewer.expand_sel_down,
            accelerator=shortcuts["Expand Selection Down"][self.os_type]
        )
        self.select_menu.add_command(
            label="Move Selected Sequences Up",
            command=self.viewer.move_sel_seqs_up,
            accelerator=shortcuts["Move Selected Sequences Up"][self.os_type]
        )
        self.select_menu.add_command(
            label="Move Selected Sequences Down",
            command=self.viewer.move_sel_seqs_down,
            accelerator=shortcuts["Move Selected Sequences Down"][self.os_type]
        )
        self.select_menu.add_command(
            label="Move Selected Sequences Top",
            command=self.viewer.move_sel_seqs_top,
            accelerator=shortcuts["Move Selected Sequences Top"][self.os_type]
        )
        self.select_menu.add_command(
            label="Move Selected Sequences Bottom",
            command=self.viewer.move_sel_seqs_bottom,
            accelerator=shortcuts["Move Selected Sequences Bottom"][self.os_type]
        )
        self.menu_bar.add_cascade(label=self.menu_names[2], menu=self.select_menu)

        # Create a "View" menu
        self.view_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.color_mode_menu = tk.Menu(self.sequence_menu, tearoff=0)
        self.color_modes = ["No Color", "Foreground Color", "Background Color"]
        self.selected_option = tk.StringVar(value=self.color_modes[0])  # default no color
        self.selected_option.trace("w", self.color_mode_selected)
        self.color_mode_menu.add_radiobutton(label=self.color_modes[0], variable=self.selected_option, value=self.color_modes[0])
        self.color_mode_menu.add_radiobutton(label=self.color_modes[1], variable=self.selected_option, value=self.color_modes[1])
        self.color_mode_menu.add_radiobutton(label=self.color_modes[2], variable=self.selected_option, value=self.color_modes[2])
        self.view_menu.add_cascade(label="Color Mode", menu=self.color_mode_menu)
        self.font_size_menu = tk.Menu(self.sequence_menu, tearoff=0)
        self.font_sizes = list(range(MIN_FONT_SIZE, MAX_FONT_SIZE + 1, 2))
        font_size = self.get_font_size()
        self.selected_font = tk.IntVar(value=font_size)  # default no color
        self.selected_font.trace("w", self.viewer.update_font_size)
        for font_size in self.font_sizes:
            self.font_size_menu.add_radiobutton(label=font_size, variable=self.selected_font, value=font_size)
        self.view_menu.add_cascade(label="Font Size", menu=self.font_size_menu)
        self.view_menu.add_command(label="Navigate to Position", command=self.viewer.navigate_to_position)
        self.menu_bar.add_cascade(label=self.menu_names[3], menu=self.view_menu)

        # Create a "Sequence" menu
        self.sequence_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.sequence_menu.add_command(
            label="Add New Sequence",
            command=self.viewer.add_new_sequence,
            accelerator=shortcuts["Add New Sequence"][self.os_type]
        )
        self.sequence_menu.add_command(
            label="Add Sequences from File",
            command=self.viewer.add_seq_from_file,
            accelerator=shortcuts["Add Sequences from File"][self.os_type]
        )
        self.sequence_menu.add_separator()
        self.sequence_menu.add_command(label="Concatenate Selected", command=self.viewer.concat_selected)
        self.comp_menu = tk.Menu(self.sequence_menu, tearoff=0)
        self.comp_menu.add_command(label="DNA", command=lambda: self.viewer.complement(rna=False))
        self.comp_menu.add_command(label="RNA", command=lambda: self.viewer.complement(rna=True))
        self.sequence_menu.add_cascade(label="Complement Selected", menu=self.comp_menu)
        self.revcomp_menu = tk.Menu(self.sequence_menu, tearoff=0)
        self.revcomp_menu.add_command(label="DNA", command=lambda: self.viewer.reverse_complement(rna=False))
        self.revcomp_menu.add_command(label="RNA", command=lambda: self.viewer.reverse_complement(rna=True))
        self.sequence_menu.add_cascade(label="Reverse Complement Selected", menu=self.revcomp_menu)
        self.transcribe_menu = tk.Menu(self.sequence_menu, tearoff=0)
        self.transcribe_menu.add_command(label="From Coding Strand",
                                         command=lambda: self.viewer.transcribe(template=False))
        self.transcribe_menu.add_command(label="From Template Strand",
                                         command=lambda: self.viewer.transcribe(template=True))
        self.sequence_menu.add_cascade(label="Transcribe Selected", menu=self.transcribe_menu)
        self.sequence_menu.add_command(label="Translate Selected", command=self.viewer.translate)
        self.sequence_menu.add_separator()
        self.menu_bar.add_cascade(label=self.menu_names[4], menu=self.sequence_menu)

        # Create a "Alignment" menu
        self.alignment_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.alignment_menu.add_command(label="Switch Alignment", command=self.switch_alignment)
        self.msa_menu = tk.Menu(self.alignment_menu, tearoff=0)
        self.msa_menu.add_command(label="Muscle", command=self.run_muscle)
        self.msa_menu.add_command(label="ClustalOmega", command=self.run_clustalomega)
        self.msa_menu.add_command(label="Prank", command=self.run_prank)
        self.msa_menu.add_command(label="Mafft", command=self.run_mafft)
        self.msa_menu.add_command(label="Probcons", command=self.run_probcons)
        self.msa_menu.add_command(label="MSAProbs", command=self.run_msaprobs)
        self.alignment_menu.add_cascade(label="Multiple Sequence Alignment", menu=self.msa_menu)
        self.pairwise_menu = tk.Menu(self.alignment_menu, tearoff=0)
        self.pairwise_menu.add_command(label="Needleman-Wunsch (Global)", command=self.run_needle)
        self.pairwise_menu.add_command(label="Smith-Waterman (Local)", command=self.run_water)
        self.alignment_menu.add_cascade(label="Pairwise Alignment", menu=self.pairwise_menu)
        self.alignment_menu.add_separator()
        self.alignment_menu.add_command(
            label="Move Selected Position Right",
            command=self.viewer.move_sel_pos_right,
            accelerator=shortcuts["Move Selected Position Right"][self.os_type]
        )
        self.alignment_menu.add_command(
            label="Move Selected Position Left",
            command=self.viewer.move_sel_pos_left,
            accelerator=shortcuts["Move Selected Position Left"][self.os_type]
        )
        self.alignment_menu.add_command(
            label="Insert Gap Move Right",
            command=self.viewer.insert_gap_move_right,
            accelerator=shortcuts["Insert Gap Move Right"][self.os_type]
        )
        self.alignment_menu.add_command(
            label="Insert Gap Move Left",
            command=self.viewer.insert_gap_move_left,
            accelerator=shortcuts["Insert Gap Move Left"][self.os_type]
        )
        self.alignment_menu.add_separator()
        self.alignment_menu.add_command(
            label="Replace Selected With Gap(s)",
            command=self.viewer.replace_sel_with_gap,
            accelerator=shortcuts["Replace Selected With Gap(s)"][self.os_type]
        )
        self.alignment_menu.add_command(
            label="Delete Gap-Only Columns",
            command=self.viewer.delete_gap_only_columns,
            accelerator=shortcuts["Delete Gap-Only Columns"][self.os_type]
        )

        # self.alignment_menu.add_command(label="Delete All Gaps", command=self.viewer.delete_all_gaps)
        self.menu_bar.add_cascade(label=self.menu_names[5], menu=self.alignment_menu)

        # Create a "Annot" menu
        self.annot_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.annot_menu.add_command(
            label="Open Selected Sequence",
            command=self.viewer.menu_open_annot,
            accelerator=shortcuts["Open Selected Sequence"][self.os_type]
        )
        self.menu_bar.add_cascade(label=self.menu_names[6], menu=self.annot_menu)

        # Create a "Tree" menu
        self.tree_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.infer_tree_menu = tk.Menu(self.tree_menu, tearoff=0)
        self.infer_tree_menu.add_command(label="FastTree", command=self.run_fasttree)
        self.infer_tree_menu.add_command(label="RAxML", command=self.run_raxml)
        self.tree_menu.add_cascade(label="Infer Tree", menu=self.infer_tree_menu)
        self.open_tree_menu = tk.Menu(self.tree_menu, tearoff=0)
        self.open_tree_menu.add_command(label="FigTree", command=self.run_fig_tree)
        self.tree_menu.add_cascade(label="Open Tree", menu=self.open_tree_menu)
        self.tree_menu.add_command(label="Convert file format", command=self.convert_tree_format)
        self.menu_bar.add_cascade(label=self.menu_names[7], menu=self.tree_menu)

        # Create a "Motif" menu
        self.motif_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.motif_menu.add_command(label="Open Motif File", command=self.read_motif)
        self.motif_menu.add_command(label="Create Motif", command=self.create_motif)
        self.menu_bar.add_cascade(label=self.menu_names[8], menu=self.motif_menu)

        # Create "Other Tools" menu
        self.other_tools_menu = tk.Menu(self.menu_bar, tearoff=0, postcommand=self.update_other_tools_menu)
        self.menu_bar.add_cascade(label=self.menu_names[9], menu=self.other_tools_menu)

    def bind_keyboard_shortcuts(self):
        # Bind for letter keys
        self.bind("<KeyPress>", self.viewer.on_key_press)

        # Bind for file_menu
        self.bind(f"<{shortcuts['New'][self.os_type]}>", lambda event: self.open_new_instance())
        self.bind(f"<{shortcuts['Open Sequence File'][self.os_type]}>", lambda event: self.open_seq_file())
        self.bind(f"<{shortcuts['Open Sequence File as Index'][self.os_type]}>",
                  lambda event: self.open_seq_file_index())
        self.bind(f"<{shortcuts['Open Alignment File'][self.os_type]}>",
                  lambda event: self.open_align_file())
        self.bind(f"<{shortcuts['Save File'][self.os_type]}>", lambda event: self.save_file())

        # Bind for edit_menu
        self.bind(f"<{shortcuts['Redo'][self.os_type]}>", lambda event: self.viewer.redo())
        self.bind(f"<{shortcuts['Copy Selection as FASTA'][self.os_type]}>",
                  lambda event: self.viewer.copy_sel_as_fasta())
        self.bind(f"<{shortcuts['Copy Selection as Text'][self.os_type]}>",
                  lambda event: self.viewer.copy_sel_as_text())
        self.bind(f"<{shortcuts['Paste FASTA'][self.os_type]}>", lambda event: self.viewer.paste_fasta())
        self.bind(f"<{shortcuts['Delete Selected'][self.os_type]}>", lambda event: self.viewer.delete_selected())
        self.bind(f"<{shortcuts['Delete Selected Letters'][self.os_type]}>",
                  lambda event: self.viewer.delete_selected_letters())
        self.bind(f"<{shortcuts['Insert Unknown Letter'][self.os_type]}>",
                  lambda event: self.viewer.insert_unknown_letter())

        # Bind for select_menu
        self.bind(f"<{shortcuts['De-select'][self.os_type]}>", lambda event: self.viewer.de_select())
        self.bind(f"<{shortcuts['Select All'][self.os_type]}>", lambda event: self.viewer.select_all())
        self.bind(f"<{shortcuts['Expand Selection Right'][self.os_type]}>",
                  lambda event: self.viewer.expand_sel_right())
        self.bind(f"<{shortcuts['Expand Selection Left'][self.os_type]}>", lambda event: self.viewer.expand_sel_left())
        self.bind(f"<{shortcuts['Expand Selection Up'][self.os_type]}>", lambda event: self.viewer.expand_sel_up())
        self.bind(f"<{shortcuts['Expand Selection Down'][self.os_type]}>", lambda event: self.viewer.expand_sel_down())
        self.bind(f"<{shortcuts['Move Selected Sequences Up'][self.os_type]}>",
                  lambda event: self.viewer.move_sel_seqs_up())
        self.bind(f"<{shortcuts['Move Selected Sequences Down'][self.os_type]}>",
                  lambda event: self.viewer.move_sel_seqs_down())
        self.bind(f"<{shortcuts['Move Selected Sequences Top'][self.os_type]}>",
                  lambda event: self.viewer.move_sel_seqs_top())
        self.bind(f"<{shortcuts['Move Selected Sequences Bottom'][self.os_type]}>",
                  lambda event: self.viewer.move_sel_seqs_bottom())

        # Bind for alignment_menu
        self.bind(f"<{shortcuts['Move Selected Position Right'][self.os_type]}>",
                  lambda event: self.viewer.move_sel_pos_right())
        self.bind(f"<{shortcuts['Move Selected Position Left'][self.os_type]}>",
                  lambda event: self.viewer.move_sel_pos_left())
        self.bind(f"<{shortcuts['Insert Gap Move Right'][self.os_type]}>",
                  lambda event: self.viewer.insert_gap_move_right())
        self.bind(f"<{shortcuts['Insert Gap Move Left'][self.os_type]}>",
                  lambda event: self.viewer.insert_gap_move_left())
        self.bind(f"<{shortcuts['Replace Selected With Gap(s)'][self.os_type]}>",
                  lambda event: self.viewer.replace_sel_with_gap())
        self.bind(f"<{shortcuts['Delete Gap-Only Columns'][self.os_type]}>",
                  lambda event: self.viewer.delete_gap_only_columns())

        # Bind for annot_menu
        self.bind(f"<{shortcuts['Open Selected Sequence'][self.os_type]}>", lambda event: self.viewer.menu_open_annot())

        # Bind for sequence_menu
        self.bind(f"<{shortcuts['Add New Sequence'][self.os_type]}>", lambda event: self.viewer.add_new_sequence())
        self.bind(f"<{shortcuts['Add Sequences from File'][self.os_type]}>",
                  lambda event: self.viewer.add_seq_from_file())

    def disable_menu_default(self):
        """
        Menu configs at start-up.( Which menu items are disabled in this mode. )
        :return:
        """

        toggle_menu_on_off(self.file_menu, ["Save File"],
                           False)

        toggle_menu_on_off(self.edit_menu, [
            "Editor Mode (Keyboard)", "Undo", "Redo", "Save Current State",
            "Copy Selection as FASTA", "Copy Selection as Text",
            "Delete Selected", "Delete Selected Letters", "Insert Unknown Letter"
        ], False)

        toggle_menu_on_off(self.select_menu, [
            "De-select", "Select All",
            "Expand Selection Right", "Expand Selection Left",
            "Expand Selection Up", "Expand Selection Down",
            "Move Selected Sequences Up", "Move Selected Sequences Down",
            "Move Selected Sequences Top", "Move Selected Sequences Bottom"
        ], False)

        toggle_menu_on_off(self.view_menu, [
            "Color Mode", "Font Size", "Navigate to Position"
        ], False)

        toggle_menu_on_off(self.sequence_menu, [
             "Concatenate Selected", "Complement Selected", "Reverse Complement Selected",
            "Transcribe Selected", "Translate Selected"
        ], False)

        toggle_menu_on_off(self.alignment_menu, [
            "Switch Alignment", "Move Selected Position Right", "Move Selected Position Left",
            "Insert Gap Move Right", "Insert Gap Move Left", "Replace Selected With Gap(s)",
            "Delete Gap-Only Columns"
        ], False)

        toggle_menu_on_off(self.annot_menu, ["Open Selected Sequence"], False)

    def disable_menu_seq(self):
        """
        Menu configs at sequence files.( Which menu items are disabled in this mode. )
        :return:
        """

        toggle_menu_on_off(self.alignment_menu, [
            "Switch Alignment", "Move Selected Position Right", "Move Selected Position Left",
            "Insert Gap Move Right", "Insert Gap Move Left", "Delete Gap-Only Columns"
        ], False)

    def disable_menu_align(self):
        """
        Menu configs at alignment files.( Which menu items are disabled in this mode. )
        :return:
        """
        toggle_menu_on_off(self.edit_menu, ["Delete Selected Letters", "Insert Unknown Letter", ], False)

        toggle_menu_on_off(self.sequence_menu, ["Concatenate Selected", "Translate Selected"], False)

    def disable_menu_seq_index(self):
        """
        Menu configs at sequence files indexing.( Which menu items are disabled in this mode. )
        :return:
        """
        toggle_menu_on_off(self.edit_menu, ["Undo", "Redo", "Paste FASTA", "Delete Selected", "Save Current State",], False)

        toggle_menu_on_off(self.select_menu, [ "Select All",
                                                "Expand Selection Right", "Expand Selection Left",
                                                "Expand Selection Up", "Expand Selection Down",
                                                "Move Selected Sequences Up", "Move Selected Sequences Down",
                                                "Move Selected Sequences Top", "Move Selected Sequences Bottom"], False)

        toggle_menu_on_off(self.view_menu, ["Font Size",], False)

        toggle_menu_on_off(self.sequence_menu, ["Add New Sequence", "Add Sequences from File",
                                                "Concatenate Selected",
                                                "Complement Selected", "Reverse Complement Selected",
                                                "Transcribe Selected", "Translate Selected",
                                                ], False)

        toggle_menu_on_off(self.alignment_menu, ["Switch Alignment",
                                                 "Move Selected Position Right",
                                                 "Move Selected Position Left",
                                                 "Insert Gap Move Right",
                                                 "Insert Gap Move Left",
                                                 "Delete Gap-Only Columns"], False)

        toggle_menu_on_off(self.motif_menu, ["Open Motif File", "Create Motif"], False)

    def disable_menu(self):
        self.menu_enabled = False

        # disable all menus
        toggle_menu_on_off_all(self.file_menu, False)
        toggle_menu_on_off_all(self.edit_menu, False)
        toggle_menu_on_off_all(self.select_menu, False)
        toggle_menu_on_off_all(self.view_menu, False)
        toggle_menu_on_off_all(self.sequence_menu, False)
        toggle_menu_on_off_all(self.alignment_menu, False)
        toggle_menu_on_off_all(self.annot_menu, False)
        toggle_menu_on_off_all(self.tree_menu, False)
        toggle_menu_on_off_all(self.motif_menu, False)
        toggle_menu_on_off_all(self.other_tools_menu, False)

        # self.other_tools_menu.configure(postcommand=None)

    def enable_menu(self):
        """
        This is not called while read thread is running. Works either before thread starts or after thread finished.
        So np for accessing self.data without lock.
        :return:
        """
        self.menu_enabled = True

        # First ENABLE ALL menus
        toggle_menu_on_off_all(self.file_menu, True)
        toggle_menu_on_off_all(self.edit_menu, True)
        toggle_menu_on_off_all(self.select_menu, True)
        toggle_menu_on_off_all(self.view_menu, True)
        toggle_menu_on_off_all(self.sequence_menu, True)
        toggle_menu_on_off_all(self.alignment_menu, True)
        toggle_menu_on_off_all(self.annot_menu, True)
        toggle_menu_on_off_all(self.tree_menu, True)
        toggle_menu_on_off_all(self.motif_menu, True)
        toggle_menu_on_off_all(self.other_tools_menu, True)

        # self.other_tools_menu.configure(postcommand=self.update_other_tools_menu)

        # Then DİSABLE BASED ON CONFİG
        if self.viewer.read_as_index:  # opened in index mode
            self.disable_menu_seq_index()
        elif self.data is not None and isinstance(self.data, Data):  # opened in seq or alignment mode
            if self.data.is_alignment:
                self.disable_menu_align()
            else:
                self.disable_menu_seq()
        else:  # no opened file
            self.disable_menu_default()

    def set_title(self, new_title):
        self.title(new_title)

    def get_font_size(self):
        font_size = self.general_settings.get("Font Size", DEFAULT_FONT_SIZE)
        font_size = int(font_size) if font_size != "" else DEFAULT_FONT_SIZE  # set to default if "" entered
        return font_size

    def check_unsaved_file(self):
        # if there is unsaved changes, refer user to save window
        save = False
        if self.unsaved_file:
            save = TKMessageBox("Warning",
                                "There are unsaved changes. Do you want to save file ?",
                                messagebox.askyesnocancel, self).show()
        return save

    def color_mode_selected(self, *args):
        if self.selected_option.get() == self.color_modes[0]:
            self.viewer.switch_color_mode(sequence_viewer.NO_COLOR)
        if self.selected_option.get() == self.color_modes[1]:
            self.viewer.switch_color_mode(sequence_viewer.FG_COLOR)
        if self.selected_option.get() == self.color_modes[2]:
            self.viewer.switch_color_mode(sequence_viewer.BG_COLOR)

        if self.viewer.read_as_index:
            self.viewer.index_controller.update_color_mode()  # paint
        else:
            self.viewer.paint_visible_lines_all()  # paint

    def open_seq_file(self):
        if self.data is not None:  # None for new instances
            self.open_new_instance(OPEN_SEQ_FILE)
            return

        # file dialog
        sfd = SequenceOpenFileDialog(self)
        if sfd.is_closed_succesfully:  # read
            file_path = sfd.in_file
            file_format = sfd.in_format
            zip_format = sfd.zip_format

            # self.print_ref_count()
            print(f"Opened File :{file_path}")

            # Create and start thread for reading
            thread = threading.Thread(target=self.read_seq_file_t, args=(file_path, file_format, zip_format))
            thread.start()

            # Popup window
            self.message_box = MessageBox(self, "Reading file...", window_title="Info", block_parent=False, ok_button=False)

            # Check read thread periodically
            self.check_thread(file_path, self.success_seq_file)

    def read_seq_file_t(self, file_path, file_format, zip_format):  # thread fun
        # time.sleep(5)
        try:
            data = SequenceReader(file_path, file_format, zip_format).read_seq_file()
            with self.read_write_lock:
                self.data = data
            print("read finished.")
        except Exception as e:
            with self.read_write_lock:
                self.data = f"Error: {e}"
            print("read error.")

    def check_thread(self, file_path, success_callback):
        # print("check_seq_file")
        # check if read thread stopped
        with self.read_write_lock:
            if self.data is not None:  # thread finished
                if isinstance(self.data, str):  # fail
                    print(self.data)
                    # inform user
                    TKMessageBox("Error", f"Could not read file. Check file format and file content.\n\n {self.data}" ,
                                 messagebox.showerror,
                                 self).show()
                    self.message_box.destroy()
                    self.data = None  # enable next open file
                else:  # success
                    success_callback(file_path)
                return  # stop scheduling

        # Schedule the next check
        self.after(100, lambda: self.check_thread(file_path, success_callback))

    def success_seq_file(self, file_path):
        if len(self.data) == 0:
            self.message_box.destroy()
            TKMessageBox("Error", f"The file does not contain any sequence data !", messagebox.showinfo, self).show()
            self.data = None  # enable next open file
        else:
            self.message_box.update_text("Displaying sequences...")
            self.update_idletasks()
            self.file_path = file_path
            self.set_title(os.path.basename(file_path))
            # time.sleep(10)
            self.viewer.display_seqs(self.data)
            self.message_box.destroy()
            # print(thread.is_alive())
            # self.print_ref_count()

    def open_seq_file_index(self):
        if self.data is not None:
            self.open_new_instance(OPEN_SEQ_FILE_INDEX)
            return

        # file dialog
        sfd = SequenceOpenIndexFileDialog(self)
        if sfd.is_closed_succesfully:  # read
            file_path = sfd.in_file
            file_format = sfd.in_format

            print(f"Opened File :{file_path}")

            # Create and start thread for reading
            thread = threading.Thread(target=self.read_seq_file_index_t, args=(file_path, file_format))
            thread.start()

            # Popup window
            self.message_box = MessageBox(self, "Reading file...", window_title="Info", block_parent=False,
                                          ok_button=False)

            # Check read thread periodically
            self.check_thread(file_path, self.success_seq_file_index)

    def read_seq_file_index_t(self, file_path, file_format):  # thread fun
        # time.sleep(5)
        try:
            cache_limit = self.general_settings.get("Index Limit", DEFAULT_INDEX_LIMIT)
            cache_limit = int(cache_limit) if cache_limit!= "" else DEFAULT_INDEX_LIMIT
            data = SequenceIndexReader(file_path, file_format, cache_limit).read_sequences()
            with self.read_write_lock:
                self.data = data
            print("read finished.")
        except Exception as e:
            with self.read_write_lock:
                self.data = f"Error: {e}"
            print("read error.")

    def success_seq_file_index(self, file_path):
        if len(self.data) == 0:
            self.message_box.destroy()
            TKMessageBox("Error", f"The file does not contain any sequence data !", messagebox.showinfo, self).show()
            self.data = None  # enable next open file
        else:
            self.message_box.update_text("Displaying sequences...")
            self.update_idletasks()
            self.file_path = file_path
            self.set_title(os.path.basename(file_path))
            # time.sleep(5)
            self.viewer.display_seqs_indexing(self.data)
            self.message_box.destroy()

    def open_align_file(self):
        if self.data is not None:
            self.open_new_instance(OPEN_ALIGN_FILE)
            return

        # file dialog
        afd = AligmentOpenFileDialog(self)
        if afd.is_closed_succesfully:  # read
            file_path = afd.in_file
            file_format = afd.in_format

            print(f"Opened File :{file_path}")

            # Create and start thread for reading
            thread = threading.Thread(target=self.read_align_file_t, args=(file_path, file_format))
            thread.start()

            # Popup window
            self.message_box = MessageBox(self, "Reading file...", window_title="Info", block_parent=False,
                                          ok_button=False)

            # Check read thread periodically
            self.check_thread(file_path, self.success_align_file)

    def read_align_file_t(self, file_path, file_format):  # thread fun
        # time.sleep(5)
        try:
            data = AlignmentReader(file_path, file_format).read_align_file()
            with self.read_write_lock:
                self.data = data
            print("read finished.")
        except Exception as e:
            with self.read_write_lock:
                self.data = f"Error: {e}"
            print("read error.")

    def success_align_file(self, file_path):
        self.message_box.destroy()
        self.file_path = file_path
        no_of_tables = self.data.number_of_tables

        if no_of_tables == 0:  # empty
            self.message_box.destroy()
            TKMessageBox("Error", f"The file does not contain any alignment data !", messagebox.showinfo, self).show()
            self.data = None   # enable next open file
            return
        elif no_of_tables == 1:
            title = os.path.basename(file_path)
        else:  # multiple - select alignment
            sa = SelectAlignment(self, no_of_tables)
            selected_table = sa.selected_table
            self.data.switch_table(selected_table - 1)
            title = os.path.basename(file_path) + f" [{selected_table}/{no_of_tables}]"

        self.message_box = MessageBox(self, "Displaying sequences...", window_title="Info", block_parent=False,
                                      ok_button=False)
        self.update_idletasks()
        self.set_title(title)
        # time.sleep(5)
        self.viewer.display_seqs(self.data)
        self.message_box.destroy()

    def save_file(self):
        if self.viewer.read_as_index:  # INDEX mode
            fd = SequenceSaveFileDialog(self)
            if fd.is_closed_succesfully:  # read
                file_path = fd.out_file
                file_format = fd.out_format
                self.save_error = None

                # Create and start thread for writing
                thread = threading.Thread(target=self.save_index_file, args=(file_path, file_format))
                thread.start()

                # Popup window
                self.message_box = MessageBox(self, "Saving file...", window_title="Info", block_parent=False,
                                              ok_button=False)

                # Check write thread periodically
                self.check_save_thread(thread)

        else:  # NORMAL mode
            if self.data.is_alignment:
                fd = AlignmentSaveFileDialog(self, self.data.number_of_tables)
            else:
                fd = SequenceSaveFileDialog(self)

            if fd.is_closed_succesfully:  # read
                file_path = fd.out_file
                file_format = fd.out_format
                sel_annots = fd.selected_annots
                self.save_error = None

                # Popup window
                self.message_box = MessageBox(self, "Saving file...", window_title="Info", block_parent=False,
                                              ok_button=False)

                # Update data according to text widget. Save current workspace.
                self.viewer.update_data()

                # Create and start thread for writing
                thread = threading.Thread(target=self.save_seq_or_align_file, args=(file_path, file_format, sel_annots))
                thread.start()

                # Check write thread periodically
                self.check_save_thread(thread)

    def save_seq_or_align_file(self, file_path, file_format, sel_annots):  # thread fun
        # time.sleep(5)
        # Write data to file with selected format.
        if self.data.is_alignment:
            writer = AlignmentWriter(self.data, file_path, file_format, sel_annots)
        else:
            writer = SequenceWriter(self.data, file_path, file_format, sel_annots)
        try:
            writer.write()
        except Exception as e:
            print(e)
            self.save_error = f"Error: {e}"
        else:
            self.save_error = None

    def save_index_file(self, file_path, file_format):  # thread fun
        # time.sleep(5)
        # Write data to file with selected format.
        writer = SequenceIndexWriter(self.data, file_path, file_format)
        try:
            writer.write()
        except Exception as e:
            print(e)
            self.save_error = f"Error: {e}"
        else:
            self.save_error = None

    def check_save_thread(self, thread):
        # print("check_save_thread")
        if not thread.is_alive():
            self.unsaved_file = False  # we saved the current file, allow user to close the app without confirming.
            self.message_box.destroy()
            if self.save_error is not None:  # thread is dead, so no lock needed
                TKMessageBox("Error", f"File couldn't be saved ! \n {self.save_error}", messagebox.showinfo, self).show()
            else:
                TKMessageBox("Info", "File saved successfully.", messagebox.showinfo, self).show()
        else:
            # Schedule the next check
            self.after(100, lambda: self.check_save_thread(thread))

    def cursor(self):
        # Get the new cursor position (index after undo)
        cursor_pos = self.viewer.seq_text_widget.index("insert")
        print(f"New cursor position: {cursor_pos}")

    def switch_alignment(self):
        no_of_tables = self.data.number_of_tables

        if no_of_tables == 1:
            return  # single table

        # select alignment
        sa = SelectAlignment(self, no_of_tables)
        selected_table = sa.selected_table - 1
        current_table = self.data.which_table()
        if current_table == selected_table:
            return  # same table

        # Update data according to text widget. Save current workspace.
        self.viewer.update_data()

        # Switch to selected table
        self.data.switch_table(selected_table)
        title = os.path.basename(self.file_path) + f" [{selected_table +1}/{no_of_tables}]"
        self.set_title(title)

        # Display new table seqs
        self.viewer.display_seqs(self.data)

        # Save current workspace at this point (after displaying new records, user cannot go back(undo/redo))
        self.viewer.save_current_state()

    def run_muscle(self):  # run in mac(v5.2: not wrapper but command runs)
        ToolWindow(self, "Alignment Tools", "Muscle")
        return

    def run_clustalomega(self):  # run in mac (v1.2.3)
        ToolWindow(self, "Alignment Tools", "ClustalOmega")
        return

    def run_prank(self):  # run in mac (v170427)
        ToolWindow(self, "Alignment Tools", "Prank")
        return

    def run_mafft(self):   # run in mac (version 7)
        ToolWindow(self, "Alignment Tools", "Mafft")
        return

    def run_probcons(self):  # run in mac (version 1.12)
        ToolWindow(self, "Alignment Tools", "Probcons")
        return

    def run_msaprobs(self):  # run on mac ( v0.9.7)
        ToolWindow(self, "Alignment Tools", "MSAProbs")
        return

    def run_needle(self):  # EMBOSS-6.6.0
        ToolWindow(self, "Alignment Tools", "Needleman-Wunsch (Global)")
        return

    def run_water(self):   # EMBOSS-6.6.0
        ToolWindow(self, "Alignment Tools", "Smith-Waterman (Local)")
        return

    def convert_align_format(self):
        caf = ConvertAlignmentFormat(self)
        if caf.is_closed_succesfully:
            in_file = caf.in_file
            in_format = caf.in_format
            out_file = caf.out_file
            out_format = caf.out_format
            self.save_error = None

            # Create and start thread for writing
            thread = threading.Thread(target=self.convert_align_format_t, args=(in_file, in_format, out_file, out_format))
            thread.start()

            # Popup window
            self.message_box = MessageBox(self, "Converting alignment file format...", window_title="Info", block_parent=False,
                                          ok_button=False)

            # Check write thread periodically
            self.check_convert_thread(thread)

    def convert_align_format_t(self, in_file, in_format, out_file, out_format):
        # time.sleep(5)
        try:
            count = AlignIO.convert(in_file, in_format, out_file, out_format)
        except Exception as e:
            print(e)
            self.save_error = f"Error: {e}"
        else:
            self.save_error = None

    def convert_seq_format(self):
        csf = ConvertSequenceFormat(self)
        if csf.is_closed_succesfully:
            in_file = csf.in_file
            in_format = csf.in_format
            out_file = csf.out_file
            out_format = csf.out_format
            self.save_error = None

            # Create and start thread for writing
            thread = threading.Thread(target=self.convert_seq_format_t, args=(in_file, in_format, out_file, out_format))
            thread.start()

            # Popup window
            self.message_box = MessageBox(self, "Converting sequence file format...", window_title="Info", block_parent=False,
                                          ok_button=False)

            # Check write thread periodically
            self.check_convert_thread(thread)

    def convert_seq_format_t(self, in_file, in_format, out_file, out_format):
        # time.sleep(5)
        try:
            count = SeqIO.convert(in_file, in_format, out_file, out_format)
        except Exception as e:
            print(e)
            self.save_error = f"Error: {e}"
        else:
            self.save_error = None

    def check_convert_thread(self, thread):
        # print("check_save_thread")
        if not thread.is_alive():
            self.message_box.destroy()
            if self.save_error is not None:  # thread is dead, so no lock needed
                TKMessageBox("Error", f"File couldn't be converted ! \n {self.save_error}", messagebox.showinfo, self).show()
            else:
                TKMessageBox("Info", "File converted successfully.", messagebox.showinfo, self).show()
        else:
            # Schedule the next check
            self.after(100, lambda: self.check_save_thread(thread))

    def open_settings(self):
        Settings(self)
        # update "other tools" menu after settings destroyed

    def run_raxml(self):  # version 8.2.12 released by Alexandros Stamatakis on May 2018
        ToolWindow(self, "Tree Tools", "RAxML")
        return

    def run_fasttree(self):  # FastTree Version 2.1.11
        ToolWindow(self, "Tree Tools", "FastTree")
        return

    def run_fig_tree(self):  # It's not working on Mac (cannot take file opens as empty), try on Windows.
        ToolWindow(self, "Tree Tools", "FigTree", os_type=self.os_type)
        return

    def convert_tree_format(self):
        ct = ConvertTree(self)
        if ct.is_closed_succesfully:
            in_file = ct.in_file
            in_format = ct.in_format
            out_file = ct.out_file
            out_format = ct.out_format
            self.save_error = None

            # Create and start thread for writing
            thread = threading.Thread(target=self.convert_tree_format_t, args=(in_file, in_format, out_file, out_format))
            thread.start()

            # Popup window
            self.message_box = MessageBox(self, "Converting tree file format...", window_title="Info", block_parent=False,
                                          ok_button=False)

            # Check write thread periodically
            self.check_convert_thread(thread)

    def convert_tree_format_t(self, in_file, in_format, out_file, out_format):
        # time.sleep(5)
        try:
            Phylo.convert(in_file, in_format, out_file, out_format)
        except Exception as e:
            print(e)
            self.save_error = f"Error: {e}"
        else:
            self.save_error = None

    def create_motif(self):
        # one dummy motif
        ms = MotifSearch(motifs=[None])

        # ask user edit mode: instances or counts
        w = SelectEditMode(self)

        if w.is_closed_succesfully:
            # get instance_mode
            is_instance_mode = w.is_instance_mode

            # show motifs
            mf = MotifFinder(self, motif_search=ms, seq_viewer=self.viewer)
            mf.create_motif(is_instance_mode)
            self.motif_finder_open = True

    def read_motif(self):
        # read from selected file
        mr = MotifFileDialog(self)
        if mr.is_closed_succesfully:
            in_file = mr.in_file
            in_format = mr.in_format
            self.motifs = None  # reset motifs

            print(f"Opened File :{in_file}")

            # Create and start thread for reading
            thread = threading.Thread(target=self.read_motif_t, args=(in_file, in_format))
            thread.start()

            # Popup window
            self.message_box = MessageBox(self, "Reading file...", window_title="Info", block_parent=False,
                                          ok_button=False)

            # Check read thread periodically
            self.check_motif(in_file)

    def read_motif_t(self, in_file, in_format):
        # time.sleep(5)
        try:
            motifs = []
            with open(in_file) as handle:
                for m in bio_motif.parse(handle, in_format):
                    motifs.append(m)
            with self.read_write_lock:
                self.motifs = motifs
            print("read finished.")
        except Exception as e:
            with self.read_write_lock:
                self.motifs = f"Error: {e}"
            print("read error.")

    def check_motif(self, in_file):
        # print("check_seq_file")
        # check if read thread stopped
        with self.read_write_lock:
            if self.motifs is not None:  # thread finished
                if isinstance(self.motifs, str):  # fail
                    print(self.motifs)
                    # inform user
                    TKMessageBox("Error", f"Could not read file. Check file format and file content.\n {self.motifs}",
                                 messagebox.showerror,
                                 self).show()
                    self.message_box.destroy()
                else:  # success
                    # time.sleep(5)
                    # create MotifSearch object and show motif finder window
                    ms = MotifSearch(motifs=self.motifs)
                    self.viewer.editor_mode_var.set(False)  # disable editor mode (replace letter will not work)
                    MotifFinder(self, motif_search=ms, seq_viewer=self.viewer)
                    self.motif_finder_open = True
                    self.message_box.destroy()

                self.motifs = None  # enable next motif (MotifSearch already got the motifs object in case of success)
                return  # stop scheduling

        # Schedule the next check
        self.after(100, lambda: self.check_motif(in_file))

    def open_new_instance(self, *args):
        """
        Two modes: bundled(pyinstaller), normal(from source code)

        STDOUT PROBLEM: (both for normal and bundled mode)

        If we do not specify the stdout argument when creating a Popen object, the new process will write its output to
        the same stdout as the parent process (the terminal or console where the Python script is executed).

        This is a problem when parent exits. When it exits, console is closed, so child cannot print to there and it is
        killed when it attempts to print anything. Therefore, we set it to DEVNULL which discards output completely.
        Note:We can open a file and set to it for later debugging.
        Note:We can apply this to the first instance (parent of all) with : sys.stdout = open("output.log", "w")

        :param args:
        :return:
        """
        if self.is_bundled:
            self.open_new_instance_bundled()  # Open Mod (args) is omitted.
        else:
            self.open_new_instance_normal(*args)

    def open_new_instance_normal(self, *args):
        # Launch a new instance
        """
        :param args:  Open Mode if created for opening new file: OPEN_SEQ_FILE, OPEN_SEQ_FILE_INDEX, OPEN_ALIGN_FILE)
        :return:
        """
        # print(sys.executable, sys.argv)

        # Launch a new instance at the calculated position (slightly right and down)
        current_x = self.winfo_x()
        current_y = self.winfo_y()
        offset_x = 30  # Horizontal offset
        offset_y = 30  # Vertical offset
        new_x = current_x + offset_x
        new_y = current_y + offset_y

        # params: exe, main script, pos, open file mode
        subprocess.Popen([sys.executable, sys.argv[0], f"--pos={new_x},{new_y}", *args],
                         stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)

    def open_new_instance_bundled(self):
        # creates empty new instance (for pyinstaller exe mode)
        # opens on top of parent
        subprocess.Popen([sys.executable], env={**os.environ, "PYINSTALLER_RESET_ENVIRONMENT": "1"},
                         stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)


if __name__ == "__main__":
    app = App('Biopyview', (600, 400))
    app.mainloop()
