from memory_profiler import profile
import threading
import tkinter as tk
from tkinter import ttk
import math

from common.configs import DEFAULT_FONT_SIZE, DEFAULT_UNDO_LIMIT
from sequence_alignment.alignment_file_dialog import AlignmentSaveFileDialog, AligmentOpenFileDialog
from sequence_alignment.alignment_reader import AlignmentReader
from sequence_editor.config import COLOR_MAPPING
from sequence_editor.sequence_file_dialog import SequenceSaveFileDialog, SequenceOpenFileDialog
from sequence_editor.sequence_reader import SequenceReader, MySeqRecord
from utils.windows import MessageBox, TKMessageBox
from viewer.annotation_viewer import AnnotationViewer
import time
from utils.model_window import ModalWindow
from utils.utils import seperator, safe_index
from utils.utils import is_valid_letter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO
import platform
from tkinter import messagebox
from io import StringIO
from common.data import Data, IndexData
import tkinter.font as tkFont

# color modes
NO_COLOR = 0
FG_COLOR = 1
BG_COLOR = 2


text_widget_configs = {
                        # "spacing1":2, "spacing2":2,  "spacing3":2,
                        "padx":1, "pady":1,
                        "bd":0, "highlightthickness":0,
                        "wrap":tk.NONE
                        }


class SequenceViewer:
    def __init__(self, master, general_settings):
        self.master = master  # window for this class' widgets (main window, App)
        self.data = None
        self.general_settings = general_settings

        self.end_time = 0
        self.start_time = 0

        # Widgets
        self.left_frame = None
        self.left_tw_frame = None
        self.seq_names_text_widget = None
        self.seq_names_ruler_widget = None
        self.right_frame = None
        self.right_tw_frame = None
        self.seq_text_widget = None
        self.seq_ruler_widget = None

        # Text selection - these control new selection
        self.text_sel_left_start_index = None   # LEFT panel
        self.text_sel_right_start_index = None  # RIGHT panel
        self.text_sel_index_start_index = None  # INDEX panel
        self.TEXT_SEL = "TS"  # used commonly by all panels

        # Status bar (this needs to be called before PanedWindow)
        self.file_type_var = None
        self.coordinate_var = None
        self.create_status_bar()

        # Create a PanedWindow and Left & Right panels
        self.paned_window = tk.PanedWindow(self.master, orient=tk.HORIZONTAL)
        self.paned_window.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.font_size = DEFAULT_FONT_SIZE  # Default font size
        self.SEQ_ID_TAG = "SEQ_ID_TAG"
        self.create_left_frame()
        self.create_right_frame()

        # Init variables
        self.editor_mode_var = tk.BooleanVar()  # for controlling only keyboard edit.
        self.edited_seq_ids = set()  # It will be added only if the sequence (Seq) changes, if its position changes, there's no need to add
        self.start_offset_in_text_widget = 1  # as text widget index counts from one.
        self.index_controller = None  # controller if file is read with indexing
        self.read_as_index = False

        # init for "add new seqs"
        self.message_box = None
        self.read_write_lock = threading.Lock()
        self.added_data = None

        # Init painting
        self.reset_prev_visible_indices()
        self.all_color_tags = [letter for letter in COLOR_MAPPING]
        self.color_mode = None
        self.last_tw_height = 0
        self.last_tw_width = 0
        self.color_off = None
        self.seqs_displayed = False
        self.switch_color_mode(NO_COLOR)  # paint initially

    def switch_color_mode(self, color_mode):
        self.color_mode = color_mode
        if color_mode == NO_COLOR:
            # delete all color tags
            # first_visible_line, last_visible_line = self.get_visible_lines()
            for tag in self.all_color_tags:
                self.seq_text_widget.tag_remove(tag, "1.0", tk.END)
            self.color_off = True
        else:
            # bg or fg
            for letter in COLOR_MAPPING:
                # tag_name = f"tag_{letter}"
                if color_mode == FG_COLOR:
                    self.seq_text_widget.tag_configure(letter, background="", foreground=COLOR_MAPPING[letter])
                if color_mode == BG_COLOR:
                    self.seq_text_widget.tag_configure(letter, background=COLOR_MAPPING[letter], foreground="")
            self.color_off = False

    def toggle_editor_mode(self):
        # print("Editor mode:", self.editor_mode_var.get())
        pass

    def create_status_bar(self):
        status_frame = tk.Frame(self.master)
        status_frame.pack(side=tk.BOTTOM, fill=tk.BOTH)
        status_frame.columnconfigure(0, weight=1, uniform="a")
        status_frame.columnconfigure(1, weight=1, uniform="a")
        status_frame.rowconfigure(0, weight=1)

        self.file_type_var = tk.StringVar(value="  File type: -")
        file_type_label = tk.Label(status_frame, font=("Arial", 12), textvariable=self.file_type_var)
        file_type_label.grid(row=0, column=0, sticky=tk.W)

        self.coordinate_var = tk.StringVar(value="Row: 0 Column: 0  ")
        coordinate_label = tk.Label(status_frame, font=("Arial", 12), textvariable=self.coordinate_var)
        coordinate_label.grid(row=0, column=1, sticky=tk.E)

    def get_font_size(self):
        font_size = self.general_settings.get("Font Size", DEFAULT_FONT_SIZE)
        font_size = int(font_size) if font_size != "" else DEFAULT_FONT_SIZE  # set to default if "" entered
        return font_size

    def get_undo_limit(self):
        undo_limit = self.general_settings.get("Undo Limit", DEFAULT_UNDO_LIMIT)
        undo_limit = int(undo_limit) if undo_limit != "" else DEFAULT_UNDO_LIMIT  # set to default if "" entered
        return undo_limit

    def create_left_frame(self):
        # Create the left frame
        self.left_frame = tk.Frame(self.paned_window)
        self.paned_window.add(self.left_frame, minsize=100)

        # Create an EMPTY Text widget
        font_size = self.get_font_size()
        self.seq_names_ruler_widget = self.create_seq_names_index_widget(self.left_frame, font_size)
        self.seq_names_ruler_widget.pack(side=tk.TOP, fill=tk.BOTH)

        # TW Frame for text widget and scrollbar ---------------------------------
        self.left_tw_frame = tk.Frame(self.left_frame)
        self.left_tw_frame.pack(side=tk.TOP, expand=True, fill=tk.BOTH)

        # Create a vertical scrollbar for the sequence names Text widget
        self.v_scrollbar_seq_names = tk.Scrollbar(self.left_tw_frame, orient=tk.VERTICAL)
        # self.v_scrollbar_seq_names.place(relx=1, rely=0, relheight=1, anchor='ne')
        self.v_scrollbar_seq_names.pack(side=tk.RIGHT, fill="y")
        self.v_scrollbar_seq_names.set(0, 1)  # this is needed for index reading and no effect for normal read.

        # Create a horizontal scrollbar for the Text widget
        self.h_scrollbar_seq_names = tk.Scrollbar(self.left_tw_frame, orient=tk.HORIZONTAL,
                                                  # command=self.seq_names_text_widget.xview
                                                  )
        # self.h_scrollbar_seq_names.place(relx=0, rely=1, relwidth=1, anchor='sw')
        self.h_scrollbar_seq_names.pack(side=tk.BOTTOM, fill="x")

        # Create a Text widget for displaying sequence names
        undo_limit = self.get_undo_limit()
        font_size = self.get_font_size()
        self.seq_names_text_widget = self.create_seq_names_text_widget(undo_limit, font_size)
        self.seq_names_text_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Scrollbar - Text widget connect
        self.h_scrollbar_seq_names.config(command=self.seq_names_text_widget.xview)
        self.seq_names_text_widget.configure(xscrollcommand=self.h_scrollbar_seq_names.set)

        # ---------------------------------------------------------------------

        # mouse configs for text widget (drag, left click etc)
        self.set_mouse_configs(self.seq_names_text_widget, self.on_mouse_drag_left_frame, self.on_mouse_release_left_frame)

        # tag for unique seq id (same color as bg not to be seen)
        self.seq_names_text_widget.tag_configure(self.SEQ_ID_TAG, foreground="white")

    def create_seq_names_index_widget(self, left_frame, font_size):
        return tk.Text(left_frame,
                       state="disabled", height=1,
                       font=("Courier", font_size),
                       **text_widget_configs)

    def create_seq_names_text_widget(self, undo_limit, font_size):
        return tk.Text(self.left_tw_frame,
                       undo=False, autoseparators=False, maxundo=undo_limit,
                       font=("Courier", font_size),
                       **text_widget_configs)

    def create_right_frame(self):
        # Create the right frame
        self.right_frame = tk.Frame(self.paned_window)
        self.paned_window.add(self.right_frame, minsize=100)  # Add to PanedWindow
        self.paned_window.sash_place(0, 200, 0)  # 200 pixels from the left edge

        # Create a Text widget for displaying index
        font_size = self.get_font_size()
        self.seq_ruler_widget = self.create_seq_index_widget(font_size)
        self.seq_ruler_widget.pack(side=tk.TOP, fill=tk.BOTH)

        # TW Frame for text widget and scrollbar -----------------------
        self.right_tw_frame = tk.Frame(self.right_frame)
        self.right_tw_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # NO VERTICAL SCROLLBAR

        # Create a horizontal scrollbar for the Text widget
        self.h_scrollbar = tk.Scrollbar(self.right_tw_frame, orient=tk.HORIZONTAL)
        # self.h_scrollbar.place(relx=0, rely=1, relwidth=1, anchor='sw')
        self.h_scrollbar.pack(side=tk.BOTTOM, fill="x")
        self.h_scrollbar.set(0, 1)  # this is needed for index reading and no effect for normal read.

        # Create a Text widget for displaying the sequence
        undo_limit = self.get_undo_limit()
        font_size = self.get_font_size()
        self.seq_text_widget = self.create_seq_text_widget(undo_limit, font_size)
        self.seq_text_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        # -----------------------------------------------------------

        # mouse configs for text widget  (drag, left click etc)
        self.set_mouse_configs(self.seq_text_widget, self.on_mouse_drag_right_frame, self.on_mouse_release_right_frame)

    def create_seq_index_widget(self, font_size):
        return tk.Text(self.right_frame,
                       height=1, state="disabled", fg="gray",
                       font=("Courier", font_size),
                       **text_widget_configs)

    def create_seq_text_widget(self, undo_limit, font_size):
        return tk.Text(self.right_tw_frame,
                       undo=False, autoseparators=False, maxundo=undo_limit,
                       font=("Courier", font_size),
                       **text_widget_configs)

    def yview_text(self, *args):
        end_time = time.time()
        runtime = end_time - self.start_time
        print("runtime", runtime)
        self.start_time = time.time()

        self.seq_names_text_widget.yview(*args)  # Update the sequence names Text widget
        self.seq_text_widget.yview(*args)  # Update the sequence Text widget

        self.display_index_line()
        self.paint_scroll_vertical()

    def xview_text(self, *args):
        self.seq_ruler_widget.xview(*args)  # Update the "index" text widget
        self.seq_text_widget.xview(*args)  # Update the "sequence" text widget

        # Note: display_index_line is not need as whole sequence is printed.
        # so horizontal scroll does not change max visible seq len
        self.paint_scroll_horizontal()

    def mouse_sync_scroll(self, event):
        # Determine platform
        sys_platform = platform.system()

        if sys_platform == 'Windows':
            # On Windows, delta is a multiple of 120
            delta = int(-1 * (event.delta / 120))
        elif sys_platform == 'Darwin':
            # On macOS, delta is a multiple of 1
            delta = int(-1 * event.delta)
        else:
            # On Linux, use <Button-4> and <Button-5>
            return

        # Scroll both text widgets when mouse wheel is used
        self.seq_names_text_widget.yview_scroll(delta, "units")
        self.seq_text_widget.yview_scroll(delta, "units")

        self.display_index_line()
        self.paint_scroll_vertical()
        return "break"

    def sync_scroll_linux(self, event):
        # On Linux, Button-4 is scroll up and Button-5 is scroll down
        if event.num == 4:
            self.seq_names_text_widget.yview_scroll(-1, "units")
            self.seq_text_widget.yview_scroll(-1, "units")
        elif event.num == 5:
            self.seq_names_text_widget.yview_scroll(1, "units")
            self.seq_text_widget.yview_scroll(1, "units")
        return "break"

    def set_mouse_wheel(self, text_widget):
        # Bind the scroll events
        text_widget.bind("<MouseWheel>", self.mouse_sync_scroll)

        # Bind for Linux systems (Button-4 and Button-5 for mouse wheel scroll)
        text_widget.bind("<Button-4>", self.sync_scroll_linux)
        text_widget.bind("<Button-5>", self.sync_scroll_linux)

    def set_mouse_configs(self, text_widget, on_mouse_drag, on_mouse_release):
        # Ignore mouse button clicks to remove default text selection/highlight.
        text_widget.bind("<Button-1>", self.ignore_mouse_event)  # Ignore left mouse button click
        text_widget.bind("<B1-Motion>", self.ignore_mouse_event)  # Ignore left mouse button drag
        text_widget.bind("<ButtonRelease-1>", self.ignore_mouse_event)  # Ignore left mouse button release

        # Square text selection/highlight.(my custom solution)
        text_widget.tag_configure(self.TEXT_SEL, background="#CDCBCA")
        text_widget.bind("<B1-Motion>", on_mouse_drag)
        text_widget.bind("<ButtonRelease-1>", on_mouse_release)

    def set_mouse_configs_null(self, text_widget):
        # Ignore mouse button clicks to remove default text selection/highlight.
        text_widget.bind("<Button-1>", self.ignore_mouse_event)  # Ignore left mouse button click
        text_widget.bind("<B1-Motion>", self.ignore_mouse_event)  # Ignore left mouse button drag
        text_widget.bind("<ButtonRelease-1>", self.ignore_mouse_event)  # Ignore left mouse button release

    def ignore_mouse_event(self, event):
        return "break"

    def on_mouse_drag_left_frame(self, event):
        # Remove selection in right panel and index panel
        self.remove_selection(self.seq_text_widget)
        self.remove_selection(self.seq_ruler_widget)

        # Get mouse cursor index. (column counting from zero and line counting from one.)
        curr_index = self.seq_names_text_widget.index("@{},{}".format(event.x, event.y))

        # Check if a "new" selection just starts: (if user just clicked for a new selection)
        # If it is not a new selection, "curr_index" is final mouse loc for previously started selection
        if self.text_sel_left_start_index is None:
            self.text_sel_left_start_index = curr_index

        # If there's a selection, remove it (as user drags, we remove old selection, add new selection)
        self.remove_selection(self.seq_names_text_widget)

        # According to final point, highlight new square region.
        start_index = list(map(int, self.text_sel_left_start_index.split('.')))
        end_index = list(map(int, curr_index.split('.')))

        first_line = min(int(start_index[0]), int(end_index[0]))
        last_line = max(int(start_index[0]), int(end_index[0]))

        for line in range(first_line, last_line + 1):
            # select
            self.select_line_left_panel(line)

            # select also right panel
            self.select_line_right_panel(line)

    def select_line_left_panel(self, line):
        highlight_len = len(self.get_seq_record_id(line))
        start_ind = f"{line}.0"
        end_ind = f"{line}.{highlight_len}"
        self.seq_names_text_widget.tag_add(self.TEXT_SEL, start_ind, end_ind)

    def select_line_right_panel(self, line):
        start_ind_right = f"{line}.0"
        end_ind_right = f"{line}.end"
        self.seq_text_widget.tag_add(self.TEXT_SEL, start_ind_right, end_ind_right)

    def on_mouse_release_left_frame(self, event):
        # This is for "left click" event (user clicked but didn't drag). We couldn't use <Button-1> event.
        self.on_mouse_drag_left_frame(event)

        # When mouse released, selection finishes.
        self.text_sel_left_start_index = None  # this starts new selection

    def on_mouse_drag_right_frame(self, event):
        # Remove selection in left and index panel
        self.remove_selection(self.seq_names_text_widget)
        self.remove_selection(self.seq_ruler_widget)

        # Get mouse cursor index. (column counting from zero and line counting from one.)
        curr_index = self.seq_text_widget.index("@{},{}".format(event.x, event.y))
        if self.read_as_index:
            self.index_controller.show_mouse_position(event)
        else:
            self.show_mouse_position(event)

        # Check if a "new" selection just starts: (if user just clicked for a new selection)
        # If it is not a new selection, "curr_index" is final mouse loc for previously started selection
        if self.text_sel_right_start_index is None:
            self.text_sel_right_start_index = curr_index

        # If there's a selection, remove it (as user drags, we remove old selection, add new selection)
        self.remove_selection(self.seq_text_widget)

        # According to final point, highlight new square region.
        start_index = list(map(int, self.text_sel_right_start_index.split('.')))
        end_index = list(map(int, curr_index.split('.')))

        first_line = min(int(start_index[0]), int(end_index[0]))
        last_line = max(int(start_index[0]), int(end_index[0]))
        first_col = min(int(start_index[1]), int(end_index[1]))
        last_col = max(int(start_index[1]), int(end_index[1]))

        for line in range(first_line, last_line + 1):
            # select
            start_ind = f"{line}.{first_col}"
            end_ind = f"{line}.{last_col + 1}"
            self.seq_text_widget.tag_add(self.TEXT_SEL, start_ind, end_ind)

            # select also left panel
            self.select_line_left_panel(line)

        self.seq_text_widget.tag_raise(self.TEXT_SEL)  # raise to top

    def on_mouse_release_right_frame(self, event):
        # This is for "left click" event (user clicked but didn't drag). We couldn't use <Button-1> event.
        self.on_mouse_drag_right_frame(event)

        # When mouse released, selection finishes.
        self.text_sel_right_start_index = None  # this starts new selection

        # self.print_selected_text()

    def on_mouse_drag_index_frame(self, event):
        # Remove selection in left and right panel
        self.remove_selection(self.seq_names_text_widget)
        self.remove_selection(self.seq_text_widget)

        # Get mouse cursor index. (column counting from zero and line counting from one.)
        curr_index = self.seq_text_widget.index("@{},{}".format(event.x, event.y))

        # Check if new selection starts: (if user just clicked for a new selection)
        # If it is not a new selection, "curr_index" is final mouse loc for previously started selection
        if self.text_sel_index_start_index is None:
            self.text_sel_index_start_index = curr_index

        # If there's a selection, remove it (as user drags, we remove old selection, add new selection)
        self.remove_selection(self.seq_ruler_widget)

        # According to final point, highlight new square region.
        start_index = list(map(int, self.text_sel_index_start_index.split('.')))
        end_index = list(map(int, curr_index.split('.')))

        # highlight index panel
        first_col = min(int(start_index[1]), int(end_index[1]))
        last_col = max(int(start_index[1]), int(end_index[1]))
        self.seq_ruler_widget.tag_add(self.TEXT_SEL, f"{1}.{first_col}", f"{1}.{last_col + 1}")

        # Get no of line in right panel
        first_line = 1
        last_line = self.get_no_of_lines(self.seq_text_widget)

        for line in range(first_line, last_line):
            start_ind = f"{line}.{first_col}"
            end_ind = f"{line}.{last_col + 1}"
            self.seq_text_widget.tag_add(self.TEXT_SEL, start_ind, end_ind)

            # select also left panel
            highlight_len = len(self.get_seq_record_id(line))
            start_ind_left = f"{line}.0"
            end_ind_left = f"{line}.{highlight_len}"
            self.seq_names_text_widget.tag_add(self.TEXT_SEL, start_ind_left, end_ind_left)

    def on_mouse_release_index_frame(self, event):
        # This is for left click event. We couldn't use <Button-1> event in the way we want.
        self.on_mouse_drag_index_frame(event)

        # When mouse released, selection finishes.
        self.text_sel_index_start_index = None  # this starts new selection

    def select_left_frame_line(self, line):
        highlight_len = len(self.get_seq_record_id(line))
        start_ind_left = f"{line}.0"
        end_ind_left = f"{line}.{highlight_len}"
        self.seq_names_text_widget.tag_add(self.TEXT_SEL, start_ind_left, end_ind_left)

    def select_right_frame_line(self, line):
        start_ind = f"{line}.0"
        end_ind = f"{line}.end"
        self.seq_text_widget.tag_add(self.TEXT_SEL, start_ind, end_ind)

    def update_font_size(self, *args):
        font_size = self.master.selected_font.get()
        increased = False
        if font_size > self.font_size:
            increased = True
        self.font_size = font_size
        self.update_font()
        if self.read_as_index:
            if increased:
                self.index_controller.font_size_increased()
            else:
                self.index_controller.font_size_decreased()
        else:
            self.display_index_line()
            self.paint_visible_lines_all()

    def update_font(self):
        """Update the font size in the text widget."""
        new_font = tkFont.Font(family="Courier", size=self.font_size)
        self.seq_names_ruler_widget.configure(font=new_font)
        self.seq_names_text_widget.configure(font=new_font)
        self.seq_ruler_widget.configure(font=new_font)
        self.seq_text_widget.configure(font=new_font)

    def remove_selection(self, text_widget):
        """
        Delete all selection in a text_widget.
        :param text_widget:
        :return:
        """
        if text_widget.tag_ranges(self.TEXT_SEL):
            text_widget.tag_remove(self.TEXT_SEL, "1.0", tk.END)

    def get_selection_bounds(self):
        """
         Each consecutive two lines are start and end of a selected line.
         This does NOT mean all lines in this region are selected.
        :return:
        """
        tag_ranges = self.seq_text_widget.tag_ranges(self.TEXT_SEL)
        if len(tag_ranges) < 2:
            return None # no selection

        pos_1 = self.seq_text_widget.index(tag_ranges[0])   # nw
        pos_2 = self.seq_text_widget.index(tag_ranges[1])   # ne
        pos_3 = self.seq_text_widget.index(tag_ranges[-2])  # sw
        pos_4 = self.seq_text_widget.index(tag_ranges[-1])  # se
        return pos_1, pos_2, pos_3, pos_4

    def get_selection(self, text_widget):
        tag_ranges = text_widget.tag_ranges(self.TEXT_SEL)
        return tag_ranges

    def get_selection_generator(self, text_widget, tag_ranges=None):
        if tag_ranges is None:
            tag_ranges = text_widget.tag_ranges(self.TEXT_SEL)

        for start, end in zip(tag_ranges[0::2], tag_ranges[1::2]):
            yield start, end

    def get_reverse_selection_generator(self, text_widget, tag_ranges=None):
        if tag_ranges is None:
            tag_ranges = text_widget.tag_ranges(self.TEXT_SEL)

        for end, start in zip(tag_ranges[::-2], tag_ranges[-2::-2]):
            yield start, end

    def any_selection(self, text_widget):
        tag_ranges = text_widget.tag_ranges(self.TEXT_SEL)
        return tag_ranges is not None and len(tag_ranges) > 0

    def print_selected_text(self):
        """
        just for debugging
        :return:
        """
        seperator("SELECTED TEXT")

        # Get all selected ranges in text widget (for re-select replaced region)
        selection = self.get_selection(self.seq_text_widget)

        # Print the start and end index of each range
        for start, end in self.get_selection_generator(self.seq_text_widget, selection):
            print(f"Start Index: {start}, End Index: {end}")

            # Get text content of the tagged range
            tagged_text = self.seq_text_widget.get(start, end)
            print("Tagged Text:", tagged_text)

            # print("all tags at start:", self.seq_text_widget.tag_names(start))

    def on_key_press(self, event):
        # Get the key symbol from the event
        # print(event.char, event.keysym)
        char = event.char
        if not self.master.motif_finder_open and is_valid_letter(char):  # ignore if motf finder open
            # replace letter
            if self.read_as_index:
                self.index_controller.replace_selected_text(char)
            else:
                self.replace_selected_text(char)
        else:
            if event.keysym in ["Left", "Right", "Up", "Down"]:
                if self.read_as_index:
                    self.index_controller.move_selection(event.keysym)
                else:
                    self.move_selection(event.keysym)

    def move_selection(self, direction):
        if not self.any_selection(self.seq_text_widget):
            return

        _, _ ,_ , start_index = self.get_selection_bounds()
        line_num, col_num = map(int, start_index.split("."))
        col_num -= 1

        # Move the tag position depending on the arrow key
        if direction == "Left":
            col_num -= 1
        elif direction == "Right":
            col_num += 1
        elif direction == "Up":
            line_num -= 1
        elif direction == "Down":
            line_num += 1

        # Construct the new position index
        new_pos = f"{line_num}.{col_num}"
        # print(new_pos)

        # Check if the new position is within valid bounds
        # if line_num > 0 and col_num >= 0 :
        letter = self.seq_text_widget.get(new_pos)
        if line_num > 0 and col_num >= 0 and letter != "\n" and letter != "":
            self.see_position(line_num, col_num, highlight_len=1)
            if direction in ["Up", "Down"]:
                self.display_index_line()
                self.paint_scroll_vertical()
            else:
                self.paint_scroll_horizontal()

    def update_edited_seqs_set(self, index):
        """
        Edited sequences are kept in a set with their unique ids.
        Given an index (line number), it finds the unique seq id and add it to the set.
        An undo may convert a seq to its original state but we cannot follow it, so if a seq is edited once it is put
        into the set.
        :param index:
        :return:
        """
        edited_line = str(index).split('.')[0]  # same for left and right panel
        print_id = self.seq_names_text_widget.get(f"{edited_line}.0", f"{edited_line}.end")
        seq_id = int(print_id.split("-")[-1])
        self.edited_seq_ids.add(seq_id)

    def check_record_annot(self, index):
        edited_line, _ = self.index_to_line_col(index)
        record, _ = self.get_line_record(edited_line)
        if isinstance(record, SeqRecord) and record.letter_annotations is not None and len(record.letter_annotations) > 0:
            return True
        else:
            return False

    def get_line_record(self, line):
        """
        given line no, returns edited info and record
        :param line:
        :return:
        """
        seq_id = self.get_seq_id(line)
        record = self.data.get_record(seq_id)
        edited = seq_id in self.edited_seq_ids
        return record, edited

    def replace_selected_text(self, replacement_char, check_editor_mode=True):
        """
        We replace only if there is selected text. So, seq_text_widget will always be edited, seq_names_text_widget not.
        :param replacement_char:
        :param check_editor_mode:
        :return:
        """
        if check_editor_mode and not self.editor_mode_var.get():
            print("Editor mode off.")
            TKMessageBox("Warning",
                         "Editor mode (Keyboard) off.",
                         messagebox.showwarning, self.master).show()
            return

        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        # Get all selected ranges in text widget (for re-select replaced region)
        selection = self.get_selection(self.seq_text_widget)

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        with DummyEditManager(self.seq_text_widget, self.seq_names_text_widget):
            # Iterate over each selected line and delete
            for start, end in self.get_selection_generator(self.seq_text_widget, tag_ranges=selection):
                # count no of letters to replace (same for each line)
                count = len(self.seq_text_widget.get(start, end))

                # Delete the characters within the tagged range
                # self.seq_text_widget.tag_remove("color", start, end)
                self.seq_text_widget.delete(start, end)

                # Insert the replacement character at the start index
                self.seq_text_widget.insert(start, count * replacement_char)

                # Add the line number to the list
                self.update_edited_seqs_set(start)

        # reselect replaced text
        for start, end in self.get_selection_generator(self.seq_text_widget, tag_ranges=selection):
            self.seq_text_widget.tag_add(self.TEXT_SEL, start, end)

        # self.seq_text_widget.update_idletasks()
        self.master.unsaved_file = True
        self.paint_visible_lines_all()

    def check_seletion_annot(self, selection):
        for start, end in self.get_selection_generator(self.seq_text_widget, tag_ranges=selection):
            if self.check_record_annot(start):
                confirm_edit = TKMessageBox("Warning",
                                            "The selection contains sequences having letter based annotations"
                                            " which will be lost after this operation. Do you want to continue ?",
                                            messagebox.askyesno, self.master).show()
                return confirm_edit
        return True

    def check_all_annot(self):
        for i in range(1, self.get_no_of_lines(self.seq_text_widget)):
            start = f"{i}.0"
            if self.check_record_annot(start):
                confirm_edit = TKMessageBox("Warning",
                                            "The file contains sequences having letter based annotations"
                                            " which will be lost after this operation.. Do you want to continue ?",
                                            messagebox.askyesno, self.master).show()
                return confirm_edit
        return True

    def navigate_to_position(self):
        n = NavigateToPosition(self.master)
        if n.is_closed_succesfully:
            try:
                # Convert to integers
                line = int(n.line)
                column = int(n.column) - 1

                # Check if line and column are greater than 1
                if not (line >= 1 and column >= 0):
                    raise ValueError
            except ValueError:
                messagebox.showerror("Invalid Input", "Line and column must be integers and greater than 1.")
                return
            else:
                if self.read_as_index:
                    self.index_controller.navigate_to_pos(line, column)
                else:
                    self.navigate_to_pos(line, column)

    def navigate_to_pos(self, line, column):
        position = self.seq_text_widget.index(f"{line}.{column}")  # bound dışındakiler bound olarak dönüyor
        self.see_position(*self.index_to_line_col(position))
        self.display_index_line()
        self.paint_visible_lines_all()

    def de_select(self):
        self.remove_selection(self.seq_text_widget)
        self.remove_selection(self.seq_names_text_widget)

    def select_all(self):
        start_time = time.time()

        number_of_lines = self.get_no_of_lines(self.seq_names_text_widget)
        for line in range(self.start_offset_in_text_widget, number_of_lines) :
            self.select_line_left_panel(line)  # select left panel
            self.select_line_right_panel(line)  # select right panel

        end_time = time.time()
        runtime = end_time - start_time
        print(f"Runtime: {runtime} seconds for select all.")

    def menu_delete_selected(self):
        self.delete_selected()

    def menu_delete_selected_letters(self):
        if self.read_as_index:
            self.index_controller.delete_selected_letters()
        else:
            self.delete_selected_letters()

    def delete_selected(self):
        """
        We delete only if there is selected "seq name". So, both tw will always be edited.
        DELETES SELECTED LINES FULLY.
        :return:
        """
        # Get all selected ranges in text widget
        selection = self.get_selection(self.seq_names_text_widget)
        if not selection:
            print("No selected text.")
            return

        with NormalEditManager(self.seq_text_widget, self.seq_names_text_widget):
            # Iterate over each tagged range in REVERSE order
            for start, end in self.get_reverse_selection_generator(self.seq_names_text_widget, selection):
                # print(start, end)
                line_no = int(str(start).split('.')[0])

                # delete fully
                self.delete_line(self.seq_text_widget, line_no)
                self.delete_line(self.seq_names_text_widget, line_no)
                # we don't add into edited list since it is deleted completely. if it comes back with undo, again np.

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def delete_selected_letters(self):
        """
        We delete only if there is selected text. So, seq_text_widget will always be edited. seq_names_text_widget
        may or may not be. We added dummy edit for seq_names_text_widget to ensure an edit is applied.
        :return:
        """
        # Get all selected ranges in text widget
        selection = self.get_selection(self.seq_text_widget)
        if not selection:
            print("No selected text.")
            return

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        with NormalEditManager(self.seq_text_widget, self.seq_names_text_widget):
            # Iterate over each tagged range in REVERSE order
            for start, end in self.get_reverse_selection_generator(self.seq_text_widget, selection):
                # print(start, end)
                line_no = int(str(start).split('.')[0])

                # check if deleted fully
                line_number, _ = self.index_to_line_col(start)
                line_len = self.get_line_len(self.seq_text_widget, line_number)
                sel_len = self.index_to_line_col(end)[1] - self.index_to_line_col(start)[1]

                if line_len == sel_len:  # deleted fully
                    self.delete_line(self.seq_text_widget, line_no)
                    self.delete_line(self.seq_names_text_widget, line_no)
                    # we don't add into edited list since it is deleted completely. if it comes back with undo, again np.
                else:
                    # Delete the characters within the tagged range
                    self.seq_text_widget.delete(start, end)

                    # Add the line number to the list
                    self.update_edited_seqs_set(start)

            # in case there is no fully deleted line, add dummy edit to seq_names_text_widget
            self.seq_names_text_widget.insert('1.0', 'x')
            self.seq_names_text_widget.delete('1.0', '1.1')

        # remove selection in seq_names_text_widget
        self.remove_selection(self.seq_names_text_widget)
        self.remove_selection(self.seq_ruler_widget)

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def check_whole_line_deleted(self, selection):
        """
        check if any of the selected lines will be fully deleted
        :param selection:
        :return:
        """
        deleted_lines = []
        for start, end in self.get_selection_generator(self.seq_text_widget, selection):
            line_number, _ = self.index_to_line_col(start)
            line_len = self.get_line_len(self.seq_text_widget, line_number)
            sel_len = int(str(end).split('.')[1]) - int(str(start).split('.')[1])
            if line_len == sel_len:
                deleted_lines.append(line_number)
        return deleted_lines

    def index_to_line_col(self, index):
        line_number = int(str(index).split('.')[0])
        col_number = int(str(index).split('.')[1])
        return line_number, col_number

    def get_line_len(self, text_widget, line_number):
        line_len = int(text_widget.index(f"{line_number}.end").split('.')[1])
        return line_len

    def insert_unknown_letter(self):
        if self.read_as_index:
            self.index_controller.insert_unknown_char()
        else:
            self.insert_unknown_char()

    def insert_unknown_char(self):
        """
          We insert only if there is selected text. So, seq_text_widget will always be edited. seq_names_text_widget not.
         :return:
         """
        if not self.any_selection(self.seq_text_widget):
            return

        selection = self.get_selection(self.seq_text_widget)

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        # add indel (?) at start of selected block
        with DummyEditManager(self.seq_text_widget, self.seq_names_text_widget):
            for start, end in self.get_selection_generator(self.seq_text_widget, selection):
                line_number = int(str(start).split('.')[0])
                self.seq_text_widget.insert(start, "?")

                # Add the line number to the list
                self.update_edited_seqs_set(start)

        self.remove_selection(self.seq_text_widget)
        self.remove_selection(self.seq_ruler_widget)
        self.remove_selection(self.seq_names_text_widget)

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def move_sel_pos_right(self):
        """
         We move only if there is indel text. So, seq_text_widget will always be edited. seq_names_text_widget not.
        :return:
        """
        if not self.any_selection(self.seq_text_widget):
            return

        selection = self.get_selection(self.seq_text_widget)

        # check if indel char "-" exist at right for all lines
        for start, end in self.get_selection_generator(self.seq_text_widget, selection):
            after_end = str(end) + "+1c"
            if self.seq_text_widget.get(end, after_end) != "-":
                return

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        # move sel right
        with DummyEditManager(self.seq_text_widget, self.seq_names_text_widget):
            for start, end in self.get_selection_generator(self.seq_text_widget, selection):
                after_start = str(start) + "+1c"
                sel_line = self.seq_text_widget.get(start, end)
                self.seq_text_widget.delete(start, end)
                self.seq_text_widget.insert(after_start, sel_line)

        # highlight moved block
        for start, end in self.get_selection_generator(self.seq_text_widget, selection):
            after_start = str(start) + "+1c"
            after_end = str(end) + "+1c"
            self.seq_text_widget.tag_add(self.TEXT_SEL, after_start, after_end)

        self.master.unsaved_file = True
        self.paint_visible_lines_all()

    def move_sel_pos_left(self):
        """
          We move only if there is indel text. So, seq_text_widget will always be edited. seq_names_text_widget not.
         :return:
         """
        if not self.any_selection(self.seq_text_widget):
            return

        selection = self.get_selection(self.seq_text_widget)

        # check if indel char "-" exist at left for all lines
        for start, end in self.get_selection_generator(self.seq_text_widget, selection):
            before_start = str(start) + "-1c"
            if self.seq_text_widget.get(before_start, start) != "-":
                return

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        # move sel left
        with DummyEditManager(self.seq_text_widget, self.seq_names_text_widget):
            for start, end in self.get_selection_generator(self.seq_text_widget, selection):
                before_start = str(start) + "-1c"
                sel_line = self.seq_text_widget.get(start, end)
                self.seq_text_widget.delete(start, end)
                self.seq_text_widget.insert(before_start, sel_line)

        # highlight moved block
        for start, end in self.get_selection_generator(self.seq_text_widget, selection):
            before_start = str(start) + "-1c"
            before_end = str(end) + "-1c"
            self.seq_text_widget.tag_add(self.TEXT_SEL, before_start, before_end)

        self.master.unsaved_file = True
        self.paint_visible_lines_all()

    def insert_gap_move_right(self):
        """
          We insert only if there is selected text. So, seq_text_widget will always be edited. seq_names_text_widget not.
         :return:
         """
        if not self.any_selection(self.seq_text_widget):
            return

        selection = self.get_selection(self.seq_text_widget)

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        # add indel (-) at start of selected block
        sel_lines = []
        with DummyEditManager(self.seq_text_widget, self.seq_names_text_widget):
            for start, end in self.get_selection_generator(self.seq_text_widget, selection):
                line_number = int(str(start).split('.')[0])
                sel_lines.append(line_number)
                self.seq_text_widget.insert(start, "-")

            # add indel at the end of line for not selected lines
            number_of_lines = self.get_no_of_lines(self.seq_text_widget)
            for line_number in range(1, number_of_lines):
                if line_number not in sel_lines:
                    self.seq_text_widget.insert(f"{line_number}.end", "-")

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def insert_gap_move_left(self):
        """
        We insert only if there is selected text. So, seq_text_widget will always be edited. seq_names_text_widget not.
        :return:
        """
        if not self.any_selection(self.seq_text_widget):
            return

        selection = self.get_selection(self.seq_text_widget)

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        # add indel (-) at start of selected block
        sel_lines = []
        with DummyEditManager(self.seq_text_widget, self.seq_names_text_widget):
            for start, end in self.get_selection_generator(self.seq_text_widget, selection):
                line_number = int(str(start).split('.')[0])
                sel_lines.append(line_number)
                self.seq_text_widget.insert(end, "-")

            # add indel at the end of line for not selected lines
            number_of_lines = self.get_no_of_lines(self.seq_text_widget)
            for line_number in range(1, number_of_lines):
                if line_number not in sel_lines:
                    self.seq_text_widget.insert(f"{line_number}.0", "-")

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def replace_sel_with_gap(self):
        if self.read_as_index:
            self.index_controller.replace_selected_text("-", check_editor_mode=False)
        else:
            self.replace_selected_text("-", check_editor_mode=False)

    def copy_sel_as_fasta(self):
        if not self.any_selection(self.seq_text_widget):
            return

        records = []
        for start, end in self.get_selection_generator(self.seq_text_widget):
            line_number = int(str(start).split('.')[0])
            seq_name_line = self.get_line(self.seq_names_text_widget, line_number)
            seq_record_id = self.parse_seq_record_id(seq_name_line)
            seq = self.seq_text_widget.get(start, end)
            seq_record = SeqRecord(seq=seq, id=seq_record_id, description="")
            records.append(seq_record)

        out_handle = StringIO()
        SeqIO.write(records, out_handle, "fasta")
        fasta_data = out_handle.getvalue()

        self.master.clipboard_clear()
        self.master.clipboard_append(fasta_data)

    def copy_sel_as_text(self):
        if not self.any_selection(self.seq_text_widget):
            return

        records = []
        for start, end in self.get_selection_generator(self.seq_text_widget):
            seq = self.seq_text_widget.get(start, end)
            records.append(seq)
            records.append("\n")

        fasta_data = "".join(records)

        self.master.clipboard_clear()
        self.master.clipboard_append(fasta_data)

    def paste_fasta(self):
        """
        If alignment is selected: add it to the current alignment, and check the length.
        :return:
        """
        if not self.check_data_availability():
            return

        # Retrieve the text from the clipboard
        clipboard_content = self.master.clipboard_get()
        new_data = self.read_fasta_from_string(clipboard_content)

        if isinstance(new_data, str):
            TKMessageBox("Error", f"Could not read file. Check file format and file content.\n {new_data}",
                         messagebox.showerror,
                         self.master).show()
            return

        if len(new_data) > 0:  # show if a file selected and read successfully
            if not self.check_alignment_size(new_data[0]):
                return

            # find visible area
            first_visible_line, last_visible_line, _, _ = self.get_visible_region()
            middle_line = (first_visible_line + last_visible_line) // 2

            self.remove_selection(self.seq_names_text_widget)
            self.remove_selection(self.seq_text_widget)

            with NormalEditManager(self.seq_names_text_widget, self.seq_text_widget):
                for seq_record in reversed(new_data):
                    new_seq_id = self.data.get_new_id()
                    self.data.add_update_record(new_seq_id, seq_record)

                    # add seq name and content to END .
                    pos = f"{middle_line}.0"
                    self.display_single_seq_name(new_seq_id, seq_record.id, pos)
                    self.display_single_seq_content(seq_record.seq, pos)

                    # select
                    self.select_left_frame_line(middle_line)
                    self.select_right_frame_line(middle_line)

            self.master.unsaved_file = True
            self.display_index_line()
            self.paint_visible_lines_all()

    def read_fasta_from_string(self, fasta_data):
        # Use StringIO to convert the string to a file-like object
        fasta_io = StringIO(fasta_data)

        try:
            # Parse the FASTA string using SeqIO or AlignIO
            if self.data.is_alignment:
                align = next(AlignIO.parse(fasta_io, "fasta"))
                sequences = list(align)
            else:
                sequences = list(SeqIO.parse(fasta_io, "fasta"))
        except Exception as e:
            sequences = f"Error: {e}"
        return sequences

    def add_seq_from_file(self):
        """
        If alignment is selected: add it to the current alignment, and check the length.
        :return:
        """
        if not self.check_data_availability():
            return

        if self.data.is_alignment:
            fd = AligmentOpenFileDialog(self.master)
        else:
            fd = SequenceOpenFileDialog(self.master)

        if fd.is_closed_succesfully:  # read
            file_path = fd.in_file
            file_format = fd.in_format
            zip_format = fd.zip_format if not self.data.is_alignment else None

            print(f"Opened File :{file_path}")

            # Create and start thread for reading
            thread = threading.Thread(target=self.read_seq_file_t, args=(file_path, file_format, zip_format))
            thread.start()

            # Popup window
            self.message_box = MessageBox(self.master, "Reading file...", window_title="Info", block_parent=False,
                                          ok_button=False)

            # Check read thread periodically
            self.check_seq_file(file_path)

    def read_seq_file_t(self, file_path, file_format, zip_format):  # thread fun
        # time.sleep(5)
        try:
            if self.data.is_alignment:
                data = AlignmentReader(file_path, file_format).read_align_file()
            else:
                data = SequenceReader(file_path, file_format, zip_format).read_seq_file()
            with self.read_write_lock:
                self.added_data = data
            print("read finished.")
        except Exception as e:
            with self.read_write_lock:
                self.added_data = f"Error: {e}"
            print("read error.")

    def check_seq_file(self, file_path):
        # print("check_seq_file")
        # check if read thread stopped
        with self.read_write_lock:
            if self.added_data is not None:  # thread finished
                if isinstance(self.added_data, str):  # fail
                    print(self.added_data)
                    # inform user
                    TKMessageBox("Error", f"Could not read file. Check file format and file content.\n\n {self.added_data}",
                                 messagebox.showerror,
                                 self.master).show()
                    self.message_box.destroy()
                    self.added_data = None  # enable next add seq
                else:  # success
                    self.message_box.update_text("Displaying sequences...")
                    self.master.update_idletasks()
                    # time.sleep(10)
                    self.display_new_seqs(self.added_data)
                    self.added_data = None  # to be collected by garbage collector + enable next add seq
                    self.message_box.destroy()
                return  # stop scheduling

        # Schedule the next check
        self.master.after(100, lambda: self.check_seq_file(file_path))

    def display_new_seqs(self, new_data):
        """
        We take the first MSA if there are multiple MSA (for alignment files)
        :param new_data:
        :return:
        """
        no_of_tables = new_data.number_of_tables
        if no_of_tables == 0 or len(new_data) == 0:  # empty file
            self.message_box.destroy()
            TKMessageBox("Error", f"The file does not contain any data !", messagebox.showinfo, self.master).show()
            return

        print("New data len: ", len(new_data))
        # check seq len
        _, first_seq_record = next(iter(new_data))
        if not self.check_alignment_size(first_seq_record):
            return

        # find visible area
        first_visible_line, last_visible_line, _, _ = self.get_visible_region()
        middle_line = (first_visible_line + last_visible_line) // 2

        self.remove_selection(self.seq_names_text_widget)
        self.remove_selection(self.seq_text_widget)

        with NormalEditManager(self.seq_names_text_widget, self.seq_text_widget):
            for _, seq_record in reversed(new_data):
                new_seq_id = self.data.get_new_id()
                self.data.add_update_record(new_seq_id, seq_record)

                # add seq name and content .
                pos = f"{middle_line}.0"
                seq_name = seq_record.id if isinstance(seq_record, SeqRecord) else seq_record.get_id()
                self.display_single_seq_name(new_seq_id, seq_name, pos)
                self.display_single_seq_content(seq_record.seq, pos)

                # select
                self.select_left_frame_line(middle_line)
                self.select_right_frame_line(middle_line)

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def add_new_sequence(self):
        if not self.check_data_availability():
            return

        w = AddNewSequence(self.master)
        if w.is_closed_succesfully:  # read
            seq_id = w.id
            seq = w.seq
            seq_record = SeqRecord(seq, id=seq_id)   # no annot.

            # check seq len
            if not self.check_alignment_size(seq_record):
                return

            # find visible area
            first_visible_line, last_visible_line, _, _ = self.get_visible_region()
            middle_line = (first_visible_line+last_visible_line)//2
            # middle_line = first_visible_line
            # middle_line = self.get_no_of_lines(self.seq_text_widget)

            with NormalEditManager(self.seq_names_text_widget, self.seq_text_widget):
                new_id = self.data.get_new_id()
                self.data.add_update_record(new_id, seq_record)

                # add seq name and content to panels.
                pos = f"{middle_line}.0"
                self.display_single_seq_name(new_id, seq_record.id, pos)
                self.display_single_seq_content(seq_record.seq, pos)

            self.remove_selection(self.seq_names_text_widget)
            self.remove_selection(self.seq_text_widget)
            self.select_left_frame_line(middle_line)
            self.select_right_frame_line(middle_line)

            self.master.unsaved_file = True
            self.display_index_line()
            self.paint_visible_lines_all()

    def check_data_availability(self):
        """
        The places where the "data" object is first created:
        - Open seq/align/index file                               --> The main creates it and passes it to the viewer.
        - Add new seq / add seq from file / paste fasta           --> The viewer creates it and passes it to the main.
        :return:
        """
        if self.data is None:  # we dont have data yet, so file type not selected.
            s = SelectFileType(self.master)
            if s.is_closed_succesfully:  # create empty data object
                self.data = Data([dict()], 0, s.is_align)
                self.master.data = self.data
                self.show_file_type(s.is_align)
            else:
                return False
        return True

    def show_file_type(self, is_align):
        file_type = "Alignment" if is_align else "Sequence"
        self.file_type_var.set(f"  File type: {file_type}")

    def show_mouse_position(self, event):
        # Get the line and column from the mouse position
        line_col = self.seq_text_widget.index(f"@{event.x},{event.y}")
        line = line_col.split(".")[0]
        col = line_col.split(".")[1]

        # Print the line and column in the format "line.column"
        # print(f"Mouse is over line.column: {line_col}")
        self.coordinate_var.set(f"Row: {line} Column: {int(col)+1}  ")

    def check_alignment_size(self, seq_record):
        """
        check seq len
        if there is no record in data, we return True
        :param seq_record:
        :return:
        """
        if self.data.is_alignment and self.get_alignment_len() != 0 and self.get_alignment_len() != len(
                seq_record):
            TKMessageBox("Warning",
                         "Sequence length is not valid. ",
                         messagebox.showerror, self.master).show()
            return False
        return True

    def get_alignment_len(self):
        # if len(self.data) > 0:
        #     seq_id, seq_record = next(iter(self.data))
        #     return len(seq_record)
        # else:
        #     return None
        return self.get_line_len(self.seq_text_widget, 1)

    def copy_seq_record(self, record:SeqRecord):
        """
        shallow copy. dont modify the internal data strs anywhere !
        That is set a new object to a field like record.seq = "..."
        Note: we change a record in two places: 1- update_data while saving file 2- annot window if save clicked.
              in both place we set fields with new objects. so shallow copy does not pose a problem.
        :param record:
        :return:
        """
        copy = SeqRecord(seq=record.seq, id=record.id, name=record.name, description=record.description,
                  dbxrefs=record.dbxrefs, features=record.features, annotations=record.annotations,
                  letter_annotations=record.letter_annotations)
        return copy

    def concat_selected(self):
        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        # Get all selected ranges in text widget
        selection = self.get_selection(self.seq_text_widget)

        # merge all selected seq records
        records = []
        for start, end in self.get_selection_generator(self.seq_text_widget, tag_ranges=selection):
            # get record
            line, col = self.index_to_line_col(start)
            record, _ = self.get_line_record(line)

            if not isinstance(record, SeqRecord):
                record = record.to_seq_record()

            # remove letter annot if not match with tw
            tw_seq = self.get_line(self.seq_text_widget, line)
            if str(record.seq) != tw_seq:
                copy = self.copy_seq_record(record)
                copy.seq = tw_seq
                copy.letter_annotations = {}
                record = copy
            records.append(record)

        if len(records) <= 1:
            print("Select at least two seqs !")
            return

        # create new seq record id for merged
        pos_1, pos_2, pos_3, pos_4 = self.get_selection_bounds()
        start_line, _ = self.index_to_line_col(pos_1)
        end_line, _ = self.index_to_line_col(pos_3)
        start_id = self.parse_seq_record_id(self.get_line(self.seq_names_text_widget, start_line))
        end_id = self.parse_seq_record_id(self.get_line(self.seq_names_text_widget, end_line))
        new_seq_record_id = start_id + "-" + end_id

        # combine records
        try:
            merged_record = records[0]  # we ensured there is at least two seqs, so merged_record cannot be same as records[0]
            for i in range(1, len(records)):
                merged_record += records[i]  # look at "SeqRecord.__add__" for annot handling
        except Exception as e:
            print(e)
            TKMessageBox("Error",
                         str(e),
                         messagebox.showerror, self.master).show()
            return

        # add merged record to data
        new_id = self.data.get_new_id()
        merged_record.id = new_seq_record_id
        self.data.add_update_record(new_id, merged_record)

        new_line = end_line + 1
        position = f"{new_line}.0"
        with NormalEditManager(self.seq_names_text_widget, self.seq_text_widget):
            self.display_single_seq_name(new_id, merged_record.id, position)
            self.display_single_seq_content(merged_record.seq, position)

        # remove selected
        self.remove_selection(self.seq_ruler_widget)
        self.remove_selection(self.seq_names_text_widget)
        self.remove_selection(self.seq_text_widget)

        # navigate to pos
        self.see_position(new_line, 0, highlight_len=len(merged_record.seq))

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def delete_gap_only_columns(self):
        """
        We delete only if there is gap only col. So, seq_text_widget will always be edited. seq_names_text_widget not.
        Only works with fixed length data (align)
        :return:
        """
        gap_columns = self.find_gap_only_columns()
        # print("gap_columns", gap_columns)
        if len(gap_columns) == 0:
            return

        # Check if letter annot exist in "any" seqs, if yes ask user to delete annot
        confirm_edit = self.check_all_annot()
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        # Sort gap-only columns in reverse order (not to affect delete indices)
        sorted_gap_columns = sorted(gap_columns, reverse=True)

        # Iterate over each line and delete gap-only columns
        first_line = self.get_line(self.seq_text_widget, 1)
        if len(gap_columns) == len(first_line):  # delete all seqs and names.
            with NormalEditManager(self.seq_text_widget, self.seq_names_text_widget):
                self.seq_names_text_widget.delete("1.0", tk.END)
                self.seq_text_widget.delete("1.0", tk.END)
        else:
            with DummyEditManager(self.seq_text_widget, self.seq_names_text_widget):
                num_lines = self.get_no_of_lines(self.seq_text_widget)
                for line_number in range(1, num_lines):
                    for column_index in sorted_gap_columns:
                        # Delete the character at the specific column index in the current line
                        start_index = f"{line_number}.{column_index}"
                        end_index = f"{line_number}.{column_index + 1}"
                        self.seq_text_widget.delete(start_index, end_index)
        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def find_gap_only_columns(self):
        # Get the total number of lines in the Text widget
        num_lines = self.get_no_of_lines(self.seq_text_widget)

        # Initialize the set with all possible indices from the first line
        first_line = self.get_line(self.seq_text_widget, 1)
        gap_columns = {i for i, char in enumerate(first_line) if char == '-'}

        # Iterate over the remaining lines
        for line_number in range(2, num_lines):
            line = self.get_line(self.seq_text_widget, line_number)
            current_gap_indices = {i for i, char in enumerate(line) if char == '-'}

            # Update the gap_columns set with intersection
            gap_columns.intersection_update(current_gap_indices)

            # If gap_columns becomes empty, break early
            if not gap_columns:
                break

        return gap_columns

    def delete_all_gaps(self):
        start_index = "1.0"
        lines_to_process = 100  # Number of lines to process in each chunk

        with DummyEditManager(self.seq_text_widget, self.seq_names_text_widget):
            while True:
                # Calculate the end index for the current chunk
                end_index = f"{int(start_index.split('.')[0]) + lines_to_process}.0"

                # Get the text from the current chunk
                chunk_text = self.seq_text_widget.get(start_index, end_index)

                if not chunk_text:
                    break  # Stop if no more text to process

                # Process the chunk
                processed_chunk = chunk_text.replace("-", "")

                # Replace the old chunk with the processed one
                self.seq_text_widget.delete(start_index, end_index)
                self.seq_text_widget.insert(start_index, processed_chunk)

                # Update the start index
                start_index = end_index

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def complement(self, rna=False):

        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        # Get all selected ranges in text widget
        selection = self.get_selection(self.seq_text_widget)
        _, _, pos_3, _ = self.get_selection_bounds()
        end_line, _ = self.index_to_line_col(pos_3)

        # get comps and store in list
        c_records = []
        for start, end in self.get_selection_generator(self.seq_text_widget, tag_ranges=selection):
            # get record
            line, _ = self.index_to_line_col(start)
            record, _ = self.get_line_record(line)

            if not isinstance(record, SeqRecord):
                record = record.to_seq_record()

            # remove letter annot if not match with tw
            tw_seq = self.get_line(self.seq_text_widget, line)
            copy = self.copy_seq_record(record)
            if str(record.seq) != tw_seq:
                copy.letter_annotations = {}
                copy.seq = Seq(tw_seq)

            # compl  (no record level complement exist)
            try:
                if rna:
                    seq = copy.seq.complement_rna()
                else:  # dna
                    seq = copy.seq.complement()
                copy.seq = seq
                comp = copy
            except Exception as e:
                print(e)
                TKMessageBox("Warning",
                             str(e),
                             messagebox.showerror, self.master).show()
                return  # exit

            # add compl record to data
            new_id = self.data.get_new_id()
            self.data.add_update_record(new_id, comp)
            c_records.append(comp)

        # determine tw pos
        new_line = end_line + 1
        pos = f"{new_line}.0"

        self.remove_selection(self.seq_names_text_widget)
        self.remove_selection(self.seq_text_widget)

        with NormalEditManager(self.seq_names_text_widget, self.seq_text_widget):
            for seq_record in reversed(c_records):
                new_seq_id = self.data.get_new_id()
                self.data.add_update_record(new_seq_id, seq_record)

                # add seq name and content .
                self.display_single_seq_name(new_seq_id, seq_record.id, pos)
                self.display_single_seq_content(seq_record.seq, pos)

                # select
                self.select_left_frame_line(new_line)
                self.select_right_frame_line(new_line)

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def reverse_complement(self, rna=False):

        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        # Get all selected ranges in text widget
        selection = self.get_selection(self.seq_text_widget)
        _, _, pos_3, _ = self.get_selection_bounds()
        end_line, _ = self.index_to_line_col(pos_3)

        # get rev comps and store in list
        rc_records = []
        for start, end in self.get_selection_generator(self.seq_text_widget, tag_ranges=selection):
            # get record
            line, _ = self.index_to_line_col(start)
            record, _ = self.get_line_record(line)

            if not isinstance(record, SeqRecord):
                record = record.to_seq_record()

            # remove letter annot if not match with tw
            tw_seq = self.get_line(self.seq_text_widget, line)
            copy = self.copy_seq_record(record)
            if str(record.seq) != tw_seq:
                copy.seq = Seq(tw_seq)
                copy.letter_annotations = {}
                copy.features = []  # it probably creates problem below (reverse_complement) if we keep it as is
            if rna:
                copy.annotations["molecule_type"] = "RNA"

            # rev compl
            try:
                # "keep the original features (with the strand and locations adjusted)."
                # "keep the original annotations (with the letter annotations reversed)."
                rev_comp = copy.reverse_complement(id=True, name=True, description=True, features=True, annotations=True,
                                          letter_annotations=True, dbxrefs=True)
            except Exception as e:
                print(e)
                TKMessageBox("Warning",
                             str(e),
                             messagebox.showerror, self.master).show()
                return  # exit

            # add rev compl record to data
            new_id = self.data.get_new_id()
            self.data.add_update_record(new_id, rev_comp)
            rc_records.append(rev_comp)

        # determine tw pos
        new_line = end_line + 1
        pos = f"{new_line}.0"

        self.remove_selection(self.seq_names_text_widget)
        self.remove_selection(self.seq_text_widget)

        with NormalEditManager(self.seq_names_text_widget, self.seq_text_widget):
            for seq_record in reversed(rc_records):
                new_seq_id = self.data.get_new_id()
                self.data.add_update_record(new_seq_id, seq_record)

                # add seq name and content .
                self.display_single_seq_name(new_seq_id, seq_record.id, pos)
                self.display_single_seq_content(seq_record.seq, pos)

                # select
                self.select_left_frame_line(new_line)
                self.select_right_frame_line(new_line)

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def transcribe(self, template=False):
        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        # Get all selected ranges in text widget
        selection = self.get_selection(self.seq_text_widget)
        _, _, pos_3, _ = self.get_selection_bounds()
        end_line, _ = self.index_to_line_col(pos_3)

        # get mRNAs and store in list
        mrna_records = []
        for start, end in self.get_selection_generator(self.seq_text_widget, tag_ranges=selection):
            # get record
            line, _ = self.index_to_line_col(start)
            record, _ = self.get_line_record(line)

            if not isinstance(record, SeqRecord):
                record = record.to_seq_record()

            # remove letter annot if not match with tw
            tw_seq = self.get_line(self.seq_text_widget, line)
            m_rna = self.copy_seq_record(record)
            if str(record.seq) != tw_seq:
                m_rna.seq = Seq(tw_seq)
                m_rna.letter_annotations = {}

            # mRNA
            try:
                if template:  # template strand (The actual biological transcription)
                    m_rna.seq = m_rna.seq.reverse_complement().transcribe()  # we used seq based reverse_complement as transcribe only avail. there

                    # Copy the old per letter annotations, "reversing" them (what about features ? )
                    letter_annotations = {}
                    for key, value in m_rna.letter_annotations.items():
                        letter_annotations[key] = value[::-1]
                    m_rna.letter_annotations = letter_annotations

                else:  # coding strand (in Biopython and bioinformatics in general,)
                    m_rna.seq = m_rna.seq.transcribe()
            except Exception as e:
                print(e)
                TKMessageBox("Warning",
                             str(e),
                             messagebox.showerror, self.master).show()
                return  # exit

            # add mRNA record to data
            new_id = self.data.get_new_id()
            self.data.add_update_record(new_id, m_rna)
            mrna_records.append(m_rna)

        # determine tw pos
        new_line = end_line + 1
        pos = f"{new_line}.0"

        self.remove_selection(self.seq_names_text_widget)
        self.remove_selection(self.seq_text_widget)

        with NormalEditManager(self.seq_names_text_widget, self.seq_text_widget):
            for seq_record in reversed(mrna_records):
                new_seq_id = self.data.get_new_id()
                self.data.add_update_record(new_seq_id, seq_record)

                # add seq name and content .
                self.display_single_seq_name(new_seq_id, seq_record.id, pos)
                self.display_single_seq_content(seq_record.seq, pos)

                # select
                self.select_left_frame_line(new_line)
                self.select_right_frame_line(new_line)

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def translate(self):
        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        # Terminate at the first in frame stop codon
        # Complete Coding Sequence (CDS)

        # codon_table = CodonTable.ambiguous_generic_by_name[table_name]
        to = TranslateOptions(self.master)
        if to.is_closed_succesfully:
            table_name = to.genetic_code
            to_stop = to.to_stop
            cds = to.cds

            # Get all selected ranges in text widget
            selection = self.get_selection(self.seq_text_widget)
            _, _, pos_3, _ = self.get_selection_bounds()
            end_line, _ = self.index_to_line_col(pos_3)

            # get aa and store in list
            aa_records = []
            for start, end in self.get_selection_generator(self.seq_text_widget, tag_ranges=selection):
                # get record
                line, _ = self.index_to_line_col(start)
                record, _ = self.get_line_record(line)

                if not isinstance(record, SeqRecord):
                    record = record.to_seq_record()

                copy = self.copy_seq_record(record)

                # update acc.to. tw
                tw_seq = self.get_line(self.seq_text_widget, line)
                if str(record.seq) != tw_seq:
                    copy.seq = tw_seq

                copy.letter_annotations = {}  # does not allow us to set True
                copy.features = []  # does not allow us to set True

                # translate
                try:
                    aa = copy.translate(table=table_name, to_stop=to_stop, cds=cds, id=True, name=True, description=True,
                                        features=False, annotations=True, letter_annotations=False, dbxrefs=True)
                except Exception as e:
                    print(e)
                    TKMessageBox("Warning",
                                 str(e),
                                 messagebox.showerror, self.master).show()
                    return  # exit

                # add aa record to data
                new_id = self.data.get_new_id()
                self.data.add_update_record(new_id, aa)
                aa_records.append(aa)

            # determine tw pos
            new_line = end_line + 1
            pos = f"{new_line}.0"

            self.remove_selection(self.seq_names_text_widget)
            self.remove_selection(self.seq_text_widget)

            with NormalEditManager(self.seq_names_text_widget, self.seq_text_widget):
                for seq_record in reversed(aa_records):
                    new_seq_id = self.data.get_new_id()
                    self.data.add_update_record(new_seq_id, seq_record)

                    # add seq name and content .
                    self.display_single_seq_name(new_seq_id, seq_record.id, pos)
                    self.display_single_seq_content(seq_record.seq, pos)

                    # select
                    self.select_left_frame_line(new_line)
                    self.select_right_frame_line(new_line)

            self.master.unsaved_file = True
            self.display_index_line()
            self.paint_visible_lines_all()

    def expand_sel_right(self):
        if not self.any_selection(self.seq_text_widget):
            return

        # Get all selected ranges in text widget (for re-select replaced region)
        selection = self.get_selection(self.seq_text_widget)

        # Print the start and end index of each range
        for start, end in self.get_selection_generator(self.seq_text_widget, selection):
            line_number, col_number = self.index_to_line_col(end)
            start_ind = f"{line_number}.{col_number}"
            end_ind = f"{line_number}.end"
            self.seq_text_widget.tag_add(self.TEXT_SEL, start_ind, end_ind)

    def expand_sel_left(self):
        if not self.any_selection(self.seq_text_widget):
            return

        # Get all selected ranges in text widget (for re-select replaced region)
        selection = self.get_selection(self.seq_text_widget)

        # Print the start and end index of each range
        for start, end in self.get_selection_generator(self.seq_text_widget, selection):
            line_number, col_number = self.index_to_line_col(start)
            start_ind = f"{line_number}.0"
            end_ind = f"{line_number}.{col_number}"
            self.seq_text_widget.tag_add(self.TEXT_SEL, start_ind, end_ind)

    def expand_sel_up(self):
        if not self.any_selection(self.seq_text_widget):
            return

        pos_1, pos_2, pos_3, pos_4 = self.get_selection_bounds()
        stop_line_number, start_col_number = self.index_to_line_col(pos_1)
        _, stop_col_number = self.index_to_line_col(pos_2)

        for line_number in range(1, stop_line_number):
            start_ind = f"{line_number}.{start_col_number}"
            end_ind = f"{line_number}.{stop_col_number}"
            self.seq_text_widget.tag_add(self.TEXT_SEL, start_ind, end_ind)

            # select also left panel
            highlight_len = len(self.get_seq_record_id(line_number))
            start_ind_left = f"{line_number}.0"
            end_ind_left = f"{line_number}.{highlight_len}"
            self.seq_names_text_widget.tag_add(self.TEXT_SEL, start_ind_left, end_ind_left)

    def expand_sel_down(self):
        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        pos_1, pos_2, pos_3, pos_4 = self.get_selection_bounds()
        start_line_number, start_col_number = self.index_to_line_col(pos_3)
        _, stop_col_number = self.index_to_line_col(pos_4)
        stop_line_number = self.get_no_of_lines(self.seq_text_widget)

        for line_number in range(start_line_number+1, stop_line_number):
            start_ind = f"{line_number}.{start_col_number}"
            end_ind = f"{line_number}.{stop_col_number}"
            self.seq_text_widget.tag_add(self.TEXT_SEL, start_ind, end_ind)

            # select also left panel
            highlight_len = len(self.get_seq_record_id(line_number))
            start_ind_left = f"{line_number}.0"
            end_ind_left = f"{line_number}.{highlight_len}"
            self.seq_names_text_widget.tag_add(self.TEXT_SEL, start_ind_left, end_ind_left)

    def move_sel_seqs_up(self):
        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        pos_1, pos_2, pos_3, pos_4 = self.get_selection_bounds()
        start_line, _ = self.index_to_line_col(pos_1)
        stop_line, _ = self.index_to_line_col(pos_3)

        if start_line == 1:
            return  # no space to go up

        with NormalEditManager(self.seq_text_widget, self.seq_names_text_widget):
            # get the line above selected region
            above_line = start_line - 1
            seq_name = self.get_line(self.seq_names_text_widget, above_line)
            seq_id = self.parse_seq_id(seq_name)
            seq_record_id = self.parse_seq_record_id(seq_name)
            seq_content = self.get_line(self.seq_text_widget, above_line)

            # put line below selected region
            self.delete_line(self.seq_names_text_widget, above_line)
            self.delete_line(self.seq_text_widget, above_line)
            self.display_single_seq_name(seq_id, seq_record_id, f"{stop_line+1}.0")
            self.display_single_seq_content(seq_content, f"{stop_line+1}.0")

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def move_sel_seqs_down(self):
        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        pos_1, pos_2, pos_3, pos_4 = self.get_selection_bounds()
        start_line, _ = self.index_to_line_col(pos_1)
        stop_line, _ = self.index_to_line_col(pos_3)

        if self.get_no_of_lines(self.seq_text_widget) - 1 == stop_line:
            return  # no space to go down

        with NormalEditManager(self.seq_text_widget, self.seq_names_text_widget):
            # get the line below selected region
            below_line = stop_line+1
            seq_name = self.get_line(self.seq_names_text_widget, below_line)
            seq_id = self.parse_seq_id(seq_name)
            seq_record_id = self.parse_seq_record_id(seq_name)
            seq_content = self.get_line(self.seq_text_widget, below_line)

            # put line above selected region
            self.delete_line(self.seq_names_text_widget, below_line)
            self.delete_line(self.seq_text_widget, below_line)
            self.display_single_seq_name(seq_id, seq_record_id, f"{start_line}.0")
            self.display_single_seq_content(seq_content, f"{start_line}.0")

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def move_sel_seqs_top(self):
        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        with NormalEditManager(self.seq_text_widget, self.seq_names_text_widget):
            # get all selected seqs
            pos_1, pos_2, pos_3, pos_4 = self.get_selection_bounds()
            start_line, _ = self.index_to_line_col(pos_1)
            stop_line, _ = self.index_to_line_col(pos_3)
            seq_names = self.seq_names_text_widget.get(f"{start_line}.0", f"{stop_line+1}.0")
            seq_contents = self.seq_text_widget.get(f"{start_line}.0", f"{stop_line+1}.0")

            # delete selected
            self.seq_names_text_widget.delete(f"{start_line}.0", f"{stop_line + 1}.0")
            self.seq_text_widget.delete(f"{start_line}.0", f"{stop_line + 1}.0")

            # insert seq names (added line by line as seq id is added with tag)
            seq_ids = [(self.parse_seq_id(seq_name), self.parse_seq_record_id(seq_name)) for seq_name in seq_names.splitlines()]
            for seq_id_pair in reversed(seq_ids):  # reverse order not to affect indices
                self.display_single_seq_name(*seq_id_pair, "1.0")

            # insert seq content
            self.seq_text_widget.insert(f"1.0", seq_contents)

            # select
            for line in range(1, stop_line-start_line+2):
                self.select_left_frame_line(line)
                self.select_right_frame_line(line)
            self.see_position(1, 0)

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def move_sel_seqs_bottom(self):
        if not self.any_selection(self.seq_text_widget):
            print("No selected text.")
            return

        with NormalEditManager(self.seq_text_widget, self.seq_names_text_widget):
            # get all selected seqs
            pos_1, pos_2, pos_3, pos_4 = self.get_selection_bounds()
            start_line, _ = self.index_to_line_col(pos_1)
            stop_line, _ = self.index_to_line_col(pos_3)
            seq_names = self.seq_names_text_widget.get(f"{start_line}.0", f"{stop_line+1}.0")
            seq_contents = self.seq_text_widget.get(f"{start_line}.0", f"{stop_line+1}.0")

            # delete selected
            self.seq_names_text_widget.delete(f"{start_line}.0", f"{stop_line + 1}.0")
            self.seq_text_widget.delete(f"{start_line}.0", f"{stop_line + 1}.0")

            bottom_index = self.seq_text_widget.index('end-1c')
            bottom_line, _ = self.index_to_line_col(bottom_index)

            # insert seq names (added line by line as seq id is added with tag)
            seq_ids = [(self.parse_seq_id(seq_name), self.parse_seq_record_id(seq_name)) for seq_name in seq_names.splitlines()]
            for seq_id_pair in reversed(seq_ids):  # reverse order not to affect indices
                self.display_single_seq_name(*seq_id_pair, bottom_index)

            # insert seq content
            self.seq_text_widget.insert(bottom_index, seq_contents)

            # select
            for line in range(bottom_line, bottom_line + (stop_line - start_line + 1)):
                self.select_left_frame_line(line)
                self.select_right_frame_line(line)
            self.see_position(bottom_line, 0)

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def save_scroll_pos(self):
        # Save the current scroll positions
        self.yview_seq = self.seq_text_widget.yview()
        self.yview_seqnames = self.seq_names_text_widget.yview()
        self.xview_seq = self.seq_text_widget.xview()
        self.xview_seqnames = self.seq_names_text_widget.xview()

    def restore_scroll_pos(self):
        # Restore the original scroll positions
        self.seq_text_widget.yview_moveto(self.yview_seq[0])
        self.seq_names_text_widget.yview_moveto(self.yview_seqnames[0])
        self.seq_text_widget.xview_moveto(self.xview_seq[0])
        self.seq_names_text_widget.xview_moveto(self.xview_seqnames[0])

    def undo(self):
        """
        Make undo for both panels.
        There is two stacks, one for undo one for redo. Edits (insert/delete) are put into "undo" stack as they are
        applied. If an "undo" is done, the edits (until reaching a seperator) are removed from "undo" stack and put into
        "redo" stack. If a "redo" is done, they are moved back to "undo" stack. If an edit is applied when there are some
        edit in "redo" stack, "redo" stack is emptied ("undo" stack is left as is) and then this edit put into "undo"
        stack, on top of existing changes. This is because we cannot make a "redo" operation when the editor is edited.

        Note: We apply undo and redo for both panels simultaneously. That is why, edits also need to progress together.
        If one of the panels is edited without other one -even once- the simultaneous behaviour ruins and cannot be
        recovered. So, if there is an edit affecting only one panel, then apply DUMMY EDIT using "make_dummy_edit" fun.

        :return:
        """
        # Enable editing
        self.seq_names_text_widget.configure(state="normal")
        self.seq_text_widget.configure(state="normal")

        # Save the current scroll positions  (tw jumps to related pos after undo-redo, we eliminate this behaviour )
        self.save_scroll_pos()

        try:
            print("tried undo")
            self.seq_names_text_widget.edit_undo()
            self.seq_text_widget.edit_undo()
        except Exception as e:
            print(f"cannot undo: {e}")
            TKMessageBox("Warning",
                         "Nothing to undo !",
                         messagebox.showwarning, self.master).show()
        else:
            print("an undo done.")
            # Restore the original scroll positions
            self.restore_scroll_pos()

            self.master.unsaved_file = True
            self.display_index_line()
            self.paint_visible_lines_all()
            self.repair_lines_without_tag()

        # Disable editing
        self.seq_names_text_widget.configure(state="disabled")
        self.seq_text_widget.configure(state="disabled")

    def redo(self):
        """
        Make redo for both panels.
        :return:
        """
        # Enable editing
        self.seq_names_text_widget.configure(state="normal")
        self.seq_text_widget.configure(state="normal")

        # Save the current scroll positions (tw jumps to related pos after undo-redo, we eliminate this behaviour )
        self.save_scroll_pos()

        try:
            print("tried redo")
            self.seq_names_text_widget.edit_redo()
            self.seq_text_widget.edit_redo()
        except Exception as e:
            print(f"cannot redo: {e}")
            TKMessageBox("Warning",
                         "Nothing to redo !",
                         messagebox.showwarning, self.master).show()
        else:
            print("a redo done.")
            # Restore the original scroll positions
            self.restore_scroll_pos()

            self.master.unsaved_file = True
            self.display_index_line()
            self.paint_visible_lines_all()
            self.repair_lines_without_tag()

        # Disable editing
        self.seq_names_text_widget.configure(state="disabled")
        self.seq_text_widget.configure(state="disabled")

    def repair_lines_without_tag(self):
        """
         after redo or undo, if a text is inserted, left panel additions do not contain seq_id tag (for white color).
         so they seem to have black color, need to repair.
        :param tag_name:
        :return:
        """
        start_time = time.time()

        total_lines = self.get_no_of_lines(self.seq_names_text_widget)

        # gather lines whose "last character" does not have the tag
        untagged_lines = [line_num for line_num in range(1, total_lines)
                          if self.SEQ_ID_TAG not in self.seq_names_text_widget.tag_names(f"{line_num}.end-1c")]

        print("untagged_lines", untagged_lines)

        for line in untagged_lines:
            line_text = self.get_line(self.seq_names_text_widget, line)
            seq_record_id = self.parse_seq_record_id(line_text)
            start_index = f"{line}.{len(seq_record_id)}"
            end_index = f"{line}.end"
            self.seq_names_text_widget.tag_add(self.SEQ_ID_TAG, start_index, end_index)

        end_time = time.time()
        runtime = end_time - start_time
        print(f"Runtime: {runtime} seconds for repair.")

    def get_line(self, text_widget, line_index):
        """
        Get line until new line char.
        :param line_index: starts from 1 ('1.0' is the position of the beginning of the text)
        :return:
        """
        line_text = text_widget.get(f"{line_index}.0", f"{line_index}.end").strip()
        return line_text

    def parse_seq_id(self, seq_name_line: str):
        """
        unique id that we defined
        :param seq_name_line:
        :return:
        """
        return int(seq_name_line.strip().split("-")[-1])

    def parse_seq_record_id(self, seq_name_line: str):
        last_dash_index = seq_name_line.rfind('-')
        return seq_name_line[:last_dash_index]

    def get_seq_record_id(self, line):
        seq_name_line = self.get_line(self.seq_names_text_widget, line)
        return self.parse_seq_record_id(seq_name_line)

    def get_seq_id(self, line):
        seq_name_line = self.get_line(self.seq_names_text_widget, line)
        return self.parse_seq_id(seq_name_line)

    def update_data(self):
        """
        Update DATA according to text widget (to pass to the writer).
        - The data will be updated based on the content in the Text widget.
        - Only the items in the "edited_seq_ids" list are checked.
        - "edited_seq_ids": Contains sequences where the Seq attribute has been modified. Sequences that have been added or deleted are not included.
        - The data will consist of the seq_ids present in the current state of the left panel.
        - In the Text widget, only the Seq attribute of a SeqRecord can be modified. Annotations (annot) can be edited in the annotation window.

        :return:
        """
        # Update common and save to disk.
        # print("Edited seq ids: ", self.edited_seq_ids)

        # # Find line indices for edited seqs. (line number corresponding to each seq id in edited seqs. )
        edited_seq_ids_l = list(self.edited_seq_ids)  # fix the order
        line_dict = self.seq_ids_to_line_number(edited_seq_ids_l)

        # Check edited seqs (if changed), update "data" accordingly. The rest is untouched.
        for seq_id, line in line_dict.items():
            # Get seq text from widget
            line_text = self.get_line(self.seq_text_widget, line)

            # Get record
            seq_record = self.data.get_record(seq_id)

            # If a seqrecord edited, create new SeqRecord object with new seq info, transfer annots.
            if str(seq_record.seq) != line_text:
                if isinstance(seq_record, SeqRecord):
                    # .letter_annotations is directly dependent to letters. So, it is not transferred.
                    new_seq_record = SeqRecord(seq=Seq(line_text), id=seq_record.id, name=seq_record.name,
                                               description=seq_record.description, annotations=seq_record.annotations,
                                               features=seq_record.features, dbxrefs=seq_record.dbxrefs)
                else:  # MySeqRecord
                    new_seq_record = MySeqRecord(seq=line_text, description=seq_record.description)

                self.data.add_update_record(seq_id, new_seq_record)  # update seq_record
                # print(f"Edited seq record ({seq_id}): ", new_seq_record)

        # reorder dict in final order of ids before saving
        final_id_order = self.get_final_id_order()
        self.data.reorder(final_id_order)

        # empty edited list
        self.edited_seq_ids.clear()

        # Save current workspace at this point
        self.save_current_state(open_popup=False)

        print("Updated data according to text widget.")

    def get_seq_ids_from_left_panel(self):  # not used -- see below fun.
        # Optimize: Read in chunks, for example, 1000 lines at a time. Split using splitlines and then parse.
        # Interaction with the GUI will be minimal, so it might be faster. Test this with a large file.
        # Compare it with fetching line by line using get.

        # Get seq ids from text widget in final order (we cannot follow order due to redo/undo operations.)
        number_of_lines = self.get_no_of_lines(self.seq_names_text_widget)
        final_id_order = [self.parse_seq_id(self.get_line(self.seq_names_text_widget, line)) for
                          line in
                          range(self.start_offset_in_text_widget, number_of_lines)]  # excludes last line, newline .
        return final_id_order

    def seq_ids_to_line_number(self, seq_ids: list):
        # optimize, like the comment above
        number_of_lines = self.get_no_of_lines(self.seq_names_text_widget)
        d = {}
        for line in range(self.start_offset_in_text_widget, number_of_lines):
            seq_id = self.parse_seq_id(self.get_line(self.seq_names_text_widget, line))
            if seq_id in seq_ids:       # get only in seq_ids
                d[seq_id] = line
            if len(d) >= len(seq_ids):  # all seq_ids found
                return d
        return d

    def get_final_id_order(self):  # returns generator
        number_of_lines = self.get_no_of_lines(self.seq_names_text_widget)
        return (self.parse_seq_id(self.get_line(self.seq_names_text_widget, line)) for
                line in range(self.start_offset_in_text_widget, number_of_lines))

    def get_no_of_lines(self, text_widget):
        """
        Includes last empty line.
        :param text_widget:
        :return:
        """
        end_index = text_widget.index('end-1c')  # this returns beginning of last empty line. for 3 seqs it returns 4.0
        number_of_lines = int(end_index.split('.')[0])
        return number_of_lines

    def save_current_state(self, open_popup=True):
        """
        Text widget is saved at this point. Cannnot go back with undo.
        :return:
        """
        #  reset undo stacks
        self.reset_undo_stacks()
        # print("Saved current state.")
        if open_popup:
            TKMessageBox("Info",
                         "Your document is saved in its current state. \n\n (Please note that you cannot revert to previous "
                         "changes using the undo feature.)",
                         messagebox.showinfo, self.master).show()

    def reset_undo_stacks(self):
        self.seq_names_text_widget.edit_reset()  # prevent undo/redo before this point.
        self.seq_text_widget.edit_reset()  # prevent undo/redo before this point.

    # @profile
    def display_seqs(self, data: Data):
        """
        Open Sequence File, Open Alignment File.
        :param data:
        :return:
        """
        self.data = data
        self.index_controller = None
        self.read_as_index = False

        # config vertical scrollbar
        self.start_time = time.time()
        self.v_scrollbar_seq_names.config(command=self.yview_text)
        self.seq_names_text_widget.configure(yscrollcommand=self.v_scrollbar_seq_names.set)

        # config horizontal scrollbar
        self.h_scrollbar.config(command=self.xview_text)
        self.seq_text_widget.configure(xscrollcommand=self.h_scrollbar.set)

        # Mouse configs for index widget (not set in index mode, only here)
        self.set_mouse_configs(self.seq_ruler_widget, self.on_mouse_drag_index_frame, self.on_mouse_release_index_frame)

        # config mouse scroll
        self.set_mouse_wheel(self.seq_names_text_widget)
        self.set_mouse_wheel(self.seq_text_widget)

        # Bind the configure event to trigger when the window size changes
        self.seq_text_widget.bind("<Configure>", self.on_resize)

        # Display the sequence names
        self.display_sequence_names()

        # Display seq contents
        self.display_all_sequences()

        # Enable undo-redo mechanism
        self.seq_text_widget.config(undo=True)
        self.seq_names_text_widget.config(undo=True)

        self.last_tw_height = self.seq_text_widget.winfo_height()

        # Display file type
        self.show_file_type(self.data.is_alignment)

        # show mouse position
        self.seq_text_widget.bind("<Motion>", self.show_mouse_position)

    def switch_undo(self):
        self.seq_text_widget.config(undo=self.seq_text_widget.cget("undo") != 1)
        self.seq_names_text_widget.config(undo=self.seq_names_text_widget.cget("undo") != 1)

    # @profile
    def display_seqs_indexing(self, data: IndexData):
        """
        Open Sequence File as Index
        :param data:
        :return:
        """
        start_time = time.time()

        self.data = data
        self.index_controller = IndexFileController(self, len(self.data), self.data)
        self.read_as_index = True

        # config vertical scrollbar
        self.v_scrollbar_seq_names.config(command=self.index_controller.yview_text)
        self.seq_names_text_widget.configure(yscrollcommand=self.index_controller.set)

        # config horizontal scrollbar
        self.h_scrollbar.config(command=self.index_controller.xview_text)
        self.seq_text_widget.configure(xscrollcommand=self.h_scrollbar.set)

        # config tw resize
        self.seq_text_widget.bind("<Configure>", self.index_controller.on_resize)

        # Mouse configs for index widget (not set in index mode, only here)
        self.set_mouse_configs_null(self.seq_ruler_widget)

        # config mouse scroll
        self.index_controller.set_mouse_wheel(self.seq_names_text_widget)
        self.index_controller.set_mouse_wheel(self.seq_text_widget)

        # display first chunk
        self.index_controller.display_sequences()

        # Display file type
        self.show_file_type(False)

        # show mouse position
        self.seq_text_widget.bind("<Motion>", self.index_controller.show_mouse_position)

        end_time = time.time()
        runtime = end_time - start_time
        print(f"Runtime: {runtime} seconds for displaying {len(self.data)} seq names.")

    def calculate_max_length(self):
        max_length = 0
        for seq_id, seq_record in self.data:
            if len(seq_record) > max_length:
                max_length = len(seq_record)
        return max_length

    def display_sequence_names(self):
        start_time = time.time()

        # Clear the text widget
        self.seq_names_text_widget.configure(state="normal")  # change to normal in order to edit text widget.
        self.seq_names_text_widget.delete("1.0", tk.END)

        # for seq_id, seq_record in self.data:
        #     if isinstance(seq_record, SeqRecord):
        #         self.display_single_seq_name(seq_id, seq_record.id, tk.END)
        #     else:  # MySeqRecord
        #         self.display_single_seq_name(seq_id, seq_record.get_id(), tk.END)

        # print seq names in chunks
        seq_name_chunk = []
        total_chunk_len = 0  # total number of letters in a chunk
        insert_threshold = 10**4  # number of letters in a chunk

        for seq_id, seq_record in self.data:
            if isinstance(seq_record, SeqRecord):
                seq_name = seq_record.id
            else:  # MySeqRecord
                seq_name = seq_record.get_id()

            seq_name_chunk.append(seq_name)
            total_chunk_len += len(seq_name)

            if total_chunk_len >= insert_threshold:
                content = "\n".join(seq_name_chunk)
                self.seq_names_text_widget.insert(tk.END, content)
                self.seq_names_text_widget.insert(tk.END, "\n")
                seq_name_chunk.clear()
                total_chunk_len = 0

        # last chunk (remaining seqs)
        if len(seq_name_chunk) > 0:
            content = "\n".join(seq_name_chunk)
            self.seq_names_text_widget.insert(tk.END, content)

            # last line
            self.seq_names_text_widget.insert(tk.END, "\n")

        # add unique ids with SEQ_ID_TAG
        for seq_id, seq_record in self.data:
            pos = self.seq_names_text_widget.index(f"{seq_id}.end")  # seq_id = line_number
            self.seq_names_text_widget.insert(pos, "-" + str(seq_id), self.SEQ_ID_TAG)

        # Update the scroll region after inserting text
        self.seq_names_text_widget.update_idletasks()

        # Disable editing.
        self.seq_names_text_widget.configure(state="disabled")

        end_time = time.time()
        runtime = end_time - start_time
        print(f"Runtime: {runtime} seconds for displaying {len(self.data)} seq names.")

    def display_single_seq_name(self, seq_id, seq_record_id, position):
        """
        Print given seq id to left panel and add new line.
        :param seq_id:  the unique integer id we defined
        :param seq_record_id:
        :param position: Should be start of a line.
        :return:
        """
        # Display seq id
        self.seq_names_text_widget.insert(position, seq_record_id)
        new_position_index = f"{position}+{len(seq_record_id)}c"   # c stands for char
        self.seq_names_text_widget.insert(new_position_index, "-" + str(seq_id) + "\n", self.SEQ_ID_TAG)

    def menu_open_annot(self):
        if self.read_as_index:
            self.index_controller.open_annot()
        else:
            self.open_annot()

    def open_annot(self):
        if self.any_selection(self.seq_text_widget):
            tw_selected = self.seq_text_widget
        elif self.any_selection(self.seq_names_text_widget):
            tw_selected = self.seq_names_text_widget
        else:
            print("no selected seq for annot")
            return

        # get first selected seq
        selection = self.get_selection(tw_selected)
        start = selection[0]  # start of first line
        line_number, _ = self.index_to_line_col(start)
        seq_id = self.get_seq_id(line_number)

        # Get seq_record
        print("Annot Seq ID :", seq_id)
        seq_record = self.data.get_record(seq_id)
        if not isinstance(seq_record, SeqRecord):  # convert MySeqRecord to SeqRecord
            seq_record = seq_record.to_seq_record()
            self.data.add_update_record(seq_id, seq_record)

        # check if seq is different in text widget
        line_text = self.get_line(self.seq_text_widget, line_number)
        show_letter_annot = (line_text == str(seq_record.seq))
        print("show_letter_annot:", seq_id, show_letter_annot)

        # we give the "seq_record" object to annot window and it updates seq_record if user click "save and close".
        # If user click "save and close", we want to save current workspace so that user cannot go back from this point
        # using undo. For that we give "save_current_workspace" method as save callback.
        AnnotationViewer(seq_id, seq_record, line_text, self.master, show_letter_annot, self.annot_seq_id_edit,
                         self.save_current_state, line_number)

    def line_to_seq_id(self, line_number):
        line= self.get_line(self.seq_names_text_widget, line_number)
        return self.parse_seq_id(line)

    def annot_seq_id_edit(self, line_number, seq_id, seq_record_id, seq_record_seq):
        """
        Runs after annot window "save and close" clicked. (no condition, run always)
        After this function is called, save_current_workspace is called in annot window and so undo stacks are emptied.
        :param line_number:
        :param seq_id:
        :param seq_record_id:
        :return:
        """
        # seq record is already updated in annot window
        if seq_id in self.edited_seq_ids:
            self.edited_seq_ids.remove(seq_id)

        with NormalEditManager(self.seq_names_text_widget, self.seq_text_widget):
            # Delete old lines
            self.delete_line(self.seq_names_text_widget, line_number)
            self.delete_line(self.seq_text_widget, line_number)

            # Insert new line
            position = f"{line_number}.0"
            self.display_single_seq_name(seq_id, seq_record_id, position)
            self.display_single_seq_content(seq_record_seq, position)

        self.master.unsaved_file = True
        self.display_index_line()
        self.paint_visible_lines_all()

    def delete_line(self, text_widget, line_number):
        """ deletes line including new line char"""
        text_widget.delete(f"{line_number}.0", f"{line_number}.0 lineend+1c")

    def add_line(self, text_widget, position, line: str):
        text_widget.insert(position, line + "\n")

    def display_index_line(self):  # Display the sequence indices line at the top
        """
        # If it has been written before, it extends or shortens it.
        # Based only on the lines in the visible region.
        :param max_length:
        :return:
        """
        max_length = self.find_max_visible_line_length()

        # Calculate the sequence indices line
        parts = []
        interval = 10
        extra_digits = 0  # extra chars due to index number instead of "." char
        start_col = 1

        # check if already exist
        index_line_len = self.get_line_len(self.seq_ruler_widget, 1)
        if index_line_len > 0:
            if max_length == index_line_len:  # index valid
                return
            elif max_length < index_line_len:  # just delete extra index
                self.seq_ruler_widget.configure(state="normal")
                self.seq_ruler_widget.delete(f"1.{max_length}", "1.end")
                self.seq_ruler_widget.configure(state="disabled")
                return
            else:  # delete last uncompleted interval and fill remaining part
                remaining = index_line_len % interval
                whole = index_line_len // interval
                chars_to_delete = remaining + 1
                self.seq_ruler_widget.configure(state="normal")
                self.seq_ruler_widget.delete(f"1.end -{chars_to_delete}c", "1.end")
                self.seq_ruler_widget.configure(state="disabled")
                start_col = max(whole * interval, 1)  # start from 1 if 0

        for i in range(start_col, max_length + 1):
            if i % interval == 0:
                index = str(i)
                remaining_space = max_length - i + 1
                index_to_write = index[:remaining_space]  # index may not fit into rem. space
                parts.append(index_to_write)
                extra_digits = len(index_to_write) - 1
            elif extra_digits > 0:   # skip this index
                extra_digits -= 1
            else:
                parts.append(".")

        seq_indices_line = "".join(parts)

        self.seq_ruler_widget.configure(state="normal")  # change to normal in order to edit text widget.
        self.seq_ruler_widget.insert(tk.END, seq_indices_line)
        self.seq_ruler_widget.configure(state="disabled")  # change to normal in order to edit text widget.

    def find_max_visible_line_length(self):
        first_visible_line, last_visible_line, = self.get_visible_lines()
        max_len = 0
        for line_number in range(first_visible_line, last_visible_line +1):
            line_len = self.get_line_len(self.seq_text_widget, line_number)
            if line_len > max_len:
                max_len = line_len
        return max_len

    def display_all_sequences(self):
        start_time = time.time()

        # Clear the Text widget
        self.seq_text_widget.configure(state="normal")  # change to normal in order to edit text widget.
        self.seq_text_widget.delete("1.0", tk.END)
        self.seq_text_widget.update_idletasks()

        # print sequences in chunks
        seq_chunk = []
        total_chunk_len = 0  # total number of letters in a chunk
        insert_threshold = 10**4  # number of letters in a chunk

        for i, (seq_id, seq_record) in enumerate(self.data):
            seq_chunk.append(str(seq_record.seq))
            total_chunk_len += len(seq_record.seq)

            if total_chunk_len >= insert_threshold:
                content = "\n".join(seq_chunk)
                self.seq_text_widget.insert(tk.END, content)
                self.seq_text_widget.insert(tk.END, "\n")
                seq_chunk.clear()
                total_chunk_len = 0

        # last chunk (remaining seqs)
        if len(seq_chunk) > 0:
            content = "\n".join(seq_chunk)
            self.seq_text_widget.insert(tk.END, content)

            # last line
            self.seq_text_widget.insert(tk.END, "\n")

        # last works
        self.display_index_line()
        self.seqs_displayed = True
        self.paint_visible_lines_all()  # Paint screen for the first time

        # Finalize
        self.seq_text_widget.update_idletasks()  # Update the scroll region after inserting text

        self.seq_text_widget.configure(state="disabled")  # Make disabled

        end_time = time.time()
        runtime = end_time - start_time
        print(f"Runtime: {runtime} seconds for displaying {len(self.data)} seqs.")

    def add_color_tag(self, seq_segment: str, line: int, first_col: int, last_col: int):
        """
        first_col and last_col included
        """
        # Loop through the sequence and group contiguous identical letters, then add corr tag
        last_letter = None
        group_start = first_col

        # while index < len(seq_segment):
        for index in range(len(seq_segment)):
            letter = seq_segment[index]
            if letter != last_letter:
                if last_letter is not None:
                    # Apply the tag for the previous group of identical letters
                    self.seq_text_widget.tag_add(last_letter, f"{line}.{group_start}", f"{line}.{first_col + index}")
                last_letter = letter
                group_start = first_col + index

        # Apply the tag for the last group of identical letters
        if last_letter is not None:
            self.seq_text_widget.tag_add(last_letter, f"{line}.{group_start}", f"{line}.{last_col + 1}")

    def paint_scroll_horizontal(self):
        """ Lines are fixed, columns are changing.. """

        # Get visible cols
        first_visible_column, last_visible_column = self.get_visible_columns()

        if not self.color_off:
            # Right or left
            # If first_visible_column or last_visible_column has increased, it has moved to the right.
            if first_visible_column > self.prev_first_visible_column or last_visible_column > self.prev_last_visible_column: # right
                uncolored_col_range = max(self.prev_last_visible_column, first_visible_column), last_visible_column
            elif first_visible_column < self.prev_first_visible_column or last_visible_column < self.prev_last_visible_column: # left
                uncolored_col_range = first_visible_column, min(self.prev_first_visible_column, last_visible_column)
            else:
                return  # If both are the same, there's no need to color. The shift is very slight, and no new letter has appeared.

            for line in range(self.prev_first_visible_line, self.prev_last_visible_line + 1):

                for tag in self.all_color_tags:
                    self.seq_text_widget.tag_remove(tag, f"{line}.{uncolored_col_range[0]}",
                                                    f"{line}.{uncolored_col_range[-1] + 1}")
                # add new color tags
                seq_segment = self.seq_text_widget.get(f"{line}.{uncolored_col_range[0]}",
                                                       f"{line}.{uncolored_col_range[-1] + 1}")
                self.add_color_tag(seq_segment, line, uncolored_col_range[0], uncolored_col_range[-1])

            # raise TEXT_SEL to top
            self.seq_text_widget.tag_raise(self.TEXT_SEL)

        # update prevs
        self.prev_first_visible_column = first_visible_column
        self.prev_last_visible_column = last_visible_column

    def paint_scroll_vertical(self):
        """ Columns are fixed, lines are changing. """
        # paint all if columns changed (which happens when horizontal scroll is at right most )
        first_visible_column, last_visible_column = self.get_visible_columns()
        if first_visible_column != self.prev_first_visible_column or last_visible_column != self.prev_last_visible_column:
            self.paint_visible_lines_all()
            return

        # paint vertical
        first_visible_line, last_visible_line = self.get_visible_lines()
        if not self.color_off:
            #up or down
            if first_visible_line > self.prev_first_visible_line or last_visible_line > self.prev_last_visible_line: # up
                uncolored_line_range = range(max(self.prev_last_visible_line, first_visible_line), last_visible_line+1)
            elif first_visible_line < self.prev_first_visible_line or last_visible_line < self.prev_last_visible_line: # down
                uncolored_line_range = range(first_visible_line, min(self.prev_first_visible_line, last_visible_line)+1)
            else:
                return  # If both are the same, there's no need to color. The shift is very slight, and no new letter has appeared.

            for line in uncolored_line_range:

                for tag in self.all_color_tags:
                    self.seq_text_widget.tag_remove(tag, f"{line}.{self.prev_first_visible_column}",
                                                    f"{line}.{self.prev_last_visible_column + 1}")

                # add new color tags
                seq_segment = self.seq_text_widget.get(f"{line}.{self.prev_first_visible_column}",
                                                       f"{line}.{self.prev_last_visible_column + 1}")
                self.add_color_tag(seq_segment, line, self.prev_first_visible_column, self.prev_last_visible_column)

            # raise TEXT_SEL to top
            self.seq_text_widget.tag_raise(self.TEXT_SEL)

        # update prevs
        self.prev_first_visible_line = first_visible_line
        self.prev_last_visible_line = last_visible_line

    def paint_resize_vertical_inc(self):
        """ Runs only if screen enlarges."""

        # paint all if columns changed (which happens when horizontal scroll is at right most )
        first_visible_column, last_visible_column = self.get_visible_columns()
        if first_visible_column != self.prev_first_visible_column or last_visible_column != self.prev_last_visible_column:
            self.paint_visible_lines_all()
            return

        # paint vertical
        first_visible_line, last_visible_line = self.get_visible_lines()
        if not self.color_off:
            # Up or down
            if first_visible_line == self.prev_first_visible_line and last_visible_line == self.prev_last_visible_line:
                return  # If both are the same, there's no need to color. The shift is very slight, and no new letter has appeared.
            elif first_visible_line == self.prev_first_visible_line:  #down
                uncolored_line_range = range(self.prev_last_visible_line, last_visible_line+1)
            elif last_visible_line == self.prev_last_visible_line:  # up
                uncolored_line_range = range(first_visible_line, self.prev_first_visible_line + 1)
            else:
                uncolored_line_range = None  # It won't enter here, this function will trigger when growing from one side.

            for line in uncolored_line_range:

                for tag in self.all_color_tags:
                    self.seq_text_widget.tag_remove(tag, f"{line}.{self.prev_first_visible_column}",
                                                    f"{line}.{self.prev_last_visible_column + 1}")

                # add new color tags
                seq_segment = self.seq_text_widget.get(f"{line}.{self.prev_first_visible_column}",
                                                       f"{line}.{self.prev_last_visible_column + 1}")
                self.add_color_tag(seq_segment, line, self.prev_first_visible_column, self.prev_last_visible_column)

            # raise TEXT_SEL to top
            self.seq_text_widget.tag_raise(self.TEXT_SEL)

        # update prevs
        self.prev_first_visible_line = first_visible_line
        self.prev_last_visible_line = last_visible_line

    def paint_resize_vertical_decr(self):
        """ Invoked only if screen shrinks."""

        # paint all if columns changed (which happens when horizontal scroll is at right most )
        first_visible_column, last_visible_column = self.get_visible_columns()
        # print(first_visible_column, last_visible_column)
        if first_visible_column != self.prev_first_visible_column or last_visible_column != self.prev_last_visible_column:
            self.paint_visible_lines_all()

    def paint_resize_horizontal(self):
        """ Invoked only if screen enlarges."""
        # Get visible cols
        first_visible_column, last_visible_column = self.get_visible_columns()

        if not self.color_off:
            # Right or left
            if first_visible_column == self.prev_first_visible_column and last_visible_column == self.prev_last_visible_column:
                return  # If both are the same, there's no need to color. The shift is very slight, and no new letter has appeared.
            elif first_visible_column == self.prev_first_visible_column:  # rigth
                uncolored_column_range = self.prev_last_visible_column, last_visible_column
            elif last_visible_column == self.prev_last_visible_column:  # left
                uncolored_column_range = first_visible_column, self.prev_first_visible_column
            else:
                uncolored_column_range = None  # It won't enter here, this function will trigger when growing from one side.

            for line in range(self.prev_first_visible_line, self.prev_last_visible_line + 1):

                for tag in self.all_color_tags:
                    self.seq_text_widget.tag_remove(tag, f"{line}.{uncolored_column_range[0]}",
                                                    f"{line}.{uncolored_column_range[-1] + 1}")

                # add new color tags
                seq_segment = self.seq_text_widget.get(f"{line}.{uncolored_column_range[0]}", f"{line}.{uncolored_column_range[-1] + 1}")
                self.add_color_tag(seq_segment, line, uncolored_column_range[0], uncolored_column_range[-1])

            # raise TEXT_SEL to top
            self.seq_text_widget.tag_raise(self.TEXT_SEL)

        # update prevs
        self.prev_first_visible_column = first_visible_column
        self.prev_last_visible_column = last_visible_column

    def update_visible_columns(self):
        """ Just updates visible columns.  """
        # Get visible cols
        first_visible_column, last_visible_column = self.get_visible_columns()

        # update prevs
        self.prev_first_visible_column = first_visible_column
        self.prev_last_visible_column = last_visible_column

    def update_visible_lines(self):
        """ Just updates visible lines.  """
        # Get visible cols
        first_visible_line, last_visible_line = self.get_visible_lines()

        # update prevs
        self.prev_first_visible_line = first_visible_line
        self.prev_last_visible_line = last_visible_line

    def update_visible_region(self):
        first_visible_line, last_visible_line = self.get_visible_lines()
        first_visible_column, last_visible_column = self.get_visible_columns()

        # update prevs
        self.prev_first_visible_line = first_visible_line
        self.prev_last_visible_line = last_visible_line

        # update prevs
        self.prev_first_visible_column = first_visible_column
        self.prev_last_visible_column = last_visible_column

    def get_visible_lines(self):
        """ Find first and last visible lines by observing SEQ_NAMES_TEXT_WİDGET. """
        first_visible_line = int(self.seq_names_text_widget.index("@0,0").split('.')[0])
        last_visible_line = int(
            self.seq_names_text_widget.index(f"@0,{self.seq_names_text_widget.winfo_height()}").split('.')[0])
        return first_visible_line, last_visible_line

    def get_visible_columns(self):  # FIRST UPDATE SEQ_INDEX_WIDGET, THEN CALL THIS METHOD !
        """ Find first and last visible lines by observing SEQ_İNDEX_WİDGET. """
        first_visible_column = int(self.seq_ruler_widget.index("@0,0").split('.')[1])
        last_visible_column = int(
            self.seq_ruler_widget.index(f"@0,{self.seq_ruler_widget.winfo_width()}").split('.')[1])
        return first_visible_column, last_visible_column

    def get_visible_region(self):  # FIRST UPDATE SEQ_INDEX_WIDGET, THEN CALL THIS METHOD !
        # Get the first and last visible line & column indices
        # Line starts from 1, while column 0.

        first_visible_line = int(self.seq_names_text_widget.index("@0,0").split('.')[0])
        last_visible_line = int(self.seq_names_text_widget.index(f"@0,{self.seq_names_text_widget.winfo_height()}").split('.')[0])
        first_visible_column = int(self.seq_ruler_widget.index("@0,0").split('.')[1])
        last_visible_column = int(self.seq_ruler_widget.index(f"@0,{self.seq_ruler_widget.winfo_width()}").split('.')[1])

        return first_visible_line, last_visible_line, first_visible_column, last_visible_column

    def paint_visible_lines_all(self):
        """ Finds visible lines and columns, paint only them.  """
        # print("paint_visible_lines_all")
        # Get visible lines and rows
        first_visible_line, last_visible_line, first_visible_column, last_visible_column = self.get_visible_region()
        # print(first_visible_line, last_visible_line, first_visible_column, last_visible_column)

        if not self.color_off:
            for line in range(first_visible_line, last_visible_line+1):
                # delete old color tags
                for tag in self.all_color_tags:
                    self.seq_text_widget.tag_remove(tag, f"{line}.{first_visible_column}",
                                                    f"{line}.{last_visible_column + 1}")

                # add new color tags
                seq_segment = self.seq_text_widget.get(f"{line}.{first_visible_column}",
                                                       f"{line}.{last_visible_column + 1}")
                self.add_color_tag(seq_segment, line, first_visible_column, last_visible_column)

            # raise TEXT_SEL to top (selected be rendered above color)
            self.seq_text_widget.tag_raise(self.TEXT_SEL)

        # update prevs
        self.prev_first_visible_line = first_visible_line
        self.prev_last_visible_line = last_visible_line
        self.prev_first_visible_column = first_visible_column
        self.prev_last_visible_column = last_visible_column

    def reset_prev_visible_indices(self):
        self.prev_first_visible_line = 0
        self.prev_last_visible_line = 0
        self.prev_first_visible_column = -1
        self.prev_last_visible_column = -1

    def on_resize(self, event):
        """
        # If the Text widget size changes, meaning the screen size changes.
        :param event:
        :return:
        """
        if not self.seqs_displayed: # Otherwise, it tries to paint during the weird resizes when the app is first created.
            return

        if event.height != self.last_tw_height:  # height changed
            # print(f"Window resized to {event.width}x{event.height}")
            self.display_index_line()

        if event.height != self.last_tw_height and event.width != self.last_tw_width:  # both changed
            self.paint_visible_lines_all()

        elif event.height > self.last_tw_height:  # heigth increased
            self.paint_resize_vertical_inc()
        elif event.height < self.last_tw_height:  # heigth decreased
            self.paint_resize_vertical_decr()

        elif event.width > self.last_tw_width:  # width increased
            self.paint_resize_horizontal()
        elif event.width < self.last_tw_width:  # width decreased -- NO PAINT, just update visible lines
            self.update_visible_columns()

        self.last_tw_height = event.height
        self.last_tw_width = event.width
        # print(f"Window resized to {event.width}x{event.height}")

    def display_single_seq_content(self, seq_record_seq, position):
        """
        Print given seq to right panel and add new line.
        :param seq_record_seq:
        :param position: Should be start of a line
        :return:
        """
        self.seq_text_widget.insert(position, str(seq_record_seq) + "\n")

    def see_position(self, line_number, col_number, highlight_len=None):

        # Create the index for the specific line
        name_start_position = f"{line_number}.0"
        seq_start_position = f"{line_number}.{col_number}"
        index_start_position = f"1.{col_number}"

        # Scroll to the specified index
        self.seq_names_text_widget.see(name_start_position)
        self.seq_text_widget.see(seq_start_position)
        self.seq_ruler_widget.see(index_start_position)

        # Highlight
        if highlight_len is not None:
            self.highlight_pos(line_number, col_number, highlight_len)

    def highlight_pos(self, line_number, col_number, highlight_len):
        name_start_position = f"{line_number}.0"
        seq_start_position = f"{line_number}.{col_number}"
        name_end_position = f"{line_number}.{len(self.get_seq_record_id(line_number))}"
        seq_end_position = f"{line_number}.{col_number + highlight_len}"

        self.remove_selection(self.seq_names_text_widget)
        self.remove_selection(self.seq_text_widget)
        self.remove_selection(self.seq_ruler_widget)

        self.seq_names_text_widget.tag_add(self.TEXT_SEL, name_start_position, name_end_position)
        self.seq_text_widget.tag_add(self.TEXT_SEL, seq_start_position, seq_end_position)

    def print_last_line(self):
        # Get the total number of lines in the Text widget
        total_lines = int(self.seq_names_text_widget.index('end-1c').split('.')[0])
        # Get the content of the last line
        for i in range(3, 0, -1):
            last_line = self.seq_names_text_widget.get(f"{total_lines-i}.0", f"{total_lines-i}.end")
            print(f"Last line: {last_line}")

def make_dummy_edit(text_widget):
    """
    dummy insertion just to trigger undo stack
    :param text_widget: of left or left panel
    :return:
    """
    text_widget.configure(state="normal")
    text_widget.edit_separator()
    text_widget.insert('1.0', 'x')
    text_widget.delete('1.0', '1.1')
    text_widget.configure(state="disabled")


class DummyEditManager:
    """
    # Note: When DummyEditManager is called, the edit_text_widget must be edited.
    # If DummyEditManager is called without editing the edit_text_widget, the two text widgets (tw) will get out of sync.
    # This is because the dummy_text_widget is always edited here, regardless.

    """
    def __init__(self, edit_text_widget, dummy_text_widget):
        """
        :param edit_text_widget:  This is edited by user.
        :param dummy_text_widget: This is edited by dummy edit.
        """
        self.edit_text_widget = edit_text_widget
        self.dummy_text_widget = dummy_text_widget

    def __enter__(self):
        # Enable editing
        self.edit_text_widget.configure(state="normal")

        # Add a separator "before" making changes (undo reverts until seperator)
        self.edit_text_widget.edit_separator()

        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        # Disable editing
        self.edit_text_widget.configure(state="disabled")

        # Call a dummy edit to left panel
        # (since this edit does not affect left panel, we need to trigger stack with dummy edit)
        make_dummy_edit(self.dummy_text_widget)


class NormalEditManager:
    """
    # Note: When NormalEditManager is called, both text widgets (tw) must be edited.
    # If one of them is not edited, the sync will be broken.
    # If neither is edited, there’s no issue (a separator is added, but there's no edit on it, and both can move until the next separator).
    # If they don't move and only the separator is deleted, the undo button will be empty, but it won't cause significant issues.

    """

    def __init__(self, first_text_widget, second_text_widget):
        self.first_text_widget = first_text_widget
        self.second_text_widget = second_text_widget

    def __enter__(self):
        # Enable editing
        self.first_text_widget.configure(state="normal")
        self.second_text_widget.configure(state="normal")

        # Add a separator "before" making changes (undo reverts until seperator)
        self.first_text_widget.edit_separator()
        self.second_text_widget.edit_separator()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        # Disable editing
        self.first_text_widget.configure(state="disabled")
        self.second_text_widget.configure(state="disabled")


class ProgressWindow(ModalWindow):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(self.parent)

        # Create window
        self.title("Progress Window")

        # self.progress_var = tk.IntVar(value=0)
        self.progress_label = tk.Label(self, text="Number of sequences: 0")
        self.progress_label.pack(pady=20)


class AddNewSequence(ModalWindow):
    def __init__(self, parent, seq_len=None):
        self.parent = parent
        super().__init__(self.parent, resizable=False)
        self.title("Add New Sequence")

        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.grid_columnconfigure(0, weight=1, uniform='a')
        self.grid_columnconfigure(1, weight=1, uniform='a')
        self.grid_columnconfigure(2, weight=1, uniform='a')
        self.grid_columnconfigure(3, weight=1, uniform='a')
        self.grid_columnconfigure(4, weight=1, uniform='a')

        self.seq_len = seq_len
        self.id = None
        self.seq = None
        self.is_closed_succesfully = False

        # Create and place labels and entry widgets
        tk.Label(self, text="ID:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.id_entry = tk.Entry(self, width=10)
        self.id_entry.grid(row=0, column=1, columnspan=4, padx=5, pady=5, sticky=tk.W)

        tk.Label(self, text="Sequence:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.seq_entry = tk.Entry(self, width=30)
        self.seq_entry.grid(row=1, column=1, columnspan=4, padx=5, pady=5, sticky=tk.W)

        # Create and place the convert button
        tk.Button(self, text="Add", command=self.add_seq).grid(row=2, column=2, padx=5, pady=20)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    def add_seq(self):
        id = self.id_entry.get().strip()
        seq = self.seq_entry.get().strip()

        if self.seq_len is not None and len(seq) != self.seq_len:
            pass
        else:
            self.id = id
            self.seq = seq
            self.is_closed_succesfully = True
            self.destroy()


class SelectFileType(ModalWindow):
    def __init__(self, parent):
        self.parent = parent
        self.is_align = False
        self.is_closed_succesfully = False
        super().__init__(self.parent, resizable=False)

        # Create a new window for editing
        self.title("Select File Type")

        # Create and place widgets
        self.label = tk.Label(self, text="Note: \"Alignment\" file type requires sequences to be fixed size.")
        self.label.pack(padx=10, pady=10)

        # self.spinbox = Spinbox(self, from_=1, to=no_of_tables)
        # self.spinbox.pack(padx=10, pady=10)

        self.opts = ["Sequence", "Alignment"]
        self.selected_var = tk.StringVar(value = self.opts[0])
        self.type_selector = tk.OptionMenu(self, self.selected_var, *self.opts)
        self.type_selector.pack(padx=10, pady=10)

        # Add a "Save" button
        self.save_button = tk.Button(self, text="Confirm", command=self.save_changes)
        self.save_button.pack(padx=10, pady=10)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    # Function to save the new values and update the Treeview
    def save_changes(self):
        self.is_align = self.selected_var.get() == self.opts[1]
        self.is_closed_succesfully = True

        self.destroy()


class NavigateToPosition(ModalWindow):
    def __init__(self, parent):
        super().__init__(parent, resizable=False)
        self.title("Navigate to Position")
        self.is_closed_succesfully = False
        self.line = None
        self.column = None

        # Create labels and entry fields for line and column
        tk.Label(self, text="Line:").grid(row=0, column=0, padx=10, pady=10)
        tk.Label(self, text="Column:").grid(row=1, column=0, padx=10, pady=10)

        self.line_entry = tk.Entry(self)
        self.column_entry = tk.Entry(self)

        self.line_entry.grid(row=0, column=1, padx=10, pady=10)
        self.column_entry.grid(row=1, column=1, padx=10, pady=10)

        # Button to navigate
        navigate_button = tk.Button(self, text="Go", command=self.navigate_to_position)
        navigate_button.grid(row=2, column=0, columnspan=2, padx=10, pady=10)

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    # Function to navigate to the specified line and column in the Text widget
    def navigate_to_position(self):
        # Get the line and column values from the entry widgets
        self.line = self.line_entry.get()
        self.column = self.column_entry.get()

        # Destroy the window
        self.is_closed_succesfully = True
        self.destroy()


class IndexFileController:
    def __init__(self, viewer: SequenceViewer, total_size, data: IndexData):
        self.viewer = viewer
        self.total_seq_size = total_size
        self.total_line_size = total_size + 1  # 1 for last empty line
        self.data = data
        self.curr_records = []  # seq records of visible lines

        self.visible_count = 0  # number of lines that fits into window, float
        self.visible_count_int = 0  # Round up visible_count. We will always print this many sequences.
        self.window_size = 0  # btw 0-1, visible_count/total_size
        self.last_f = 0  # last position btw 0-1
        self.first_line_int = 1  # The first sequence printed in the tw is actually which sequence (based on total_seq_size).
        self.last_line_int = 1  # The last sequence printed in the tw is actually which sequence (based on total_seq_size).
        # self.last_offset = 0  # How much the first visible line is offset.
        self.hit_bottom = False
        self.hit_top = False
        self.stop_deleting = False
        self.last_tw_height = 0
        self.last_tw_width = 0

        # paint visible area
        self.prev_first_visible_column = -1
        self.prev_last_visible_column = -1

    def display_sequences(self):
        """
        Only runs at START
        :return:
        """
        self.update_state()
        print("visible_count", self.visible_count)
        self.moveto(0)  # by default it is in NO COLOR mode, so it won't attempt to paint.
        # initialize visible columns (we need index line. it is filled in moveto )
        self.prev_first_visible_column, self.prev_last_visible_column = self.viewer.get_visible_columns()

        self.viewer.seqs_displayed = True
        self.last_tw_height = self.viewer.seq_text_widget.winfo_height()

    def update_state(self):
        self.update_visible_count()
        self.update_window_size()

    def update_visible_count(self):
        """
        Returns the exact number of lines that fit into the visible area as a float.
        """

        # Get the height of the first line in pixels (we used index line as it remains always and they are same)
        ruler_line_metrics = self.viewer.seq_ruler_widget.dlineinfo('1.0')  # (padx, pady, width, height, ?) # It will only change when the font size changes. It can be optimized.
        widget_height = self.viewer.seq_names_text_widget.winfo_height() - ruler_line_metrics[1]   # subtract pady (space at top)

        self.visible_count = 1
        if ruler_line_metrics:
            line_height = ruler_line_metrics[3]
            if line_height > 0:
                self.visible_count = widget_height / line_height  # Float value representing exact number of lines

        self.visible_count_int = math.ceil(self.visible_count) + 1  # It would be one more than the visible count. Maybe we don’t need +1 since the offset is removed.

        # print("visible_count:", self.visible_count)
        # print("visible lines:", self.viewer.get_visible_lines())

    def update_window_size(self):
        self.window_size = self.visible_count / self.total_line_size

    def xview_text(self, *args):
        self.viewer.seq_ruler_widget.xview(*args)  # Update the "index" text widget
        self.viewer.seq_text_widget.xview(*args)  # Update the "sequence" text widget

        # Note: display_index_line is not need as whole sequence is printed.
        # so horizontal scroll does not change max visible seq len
        self.paint_scroll_horizontal()

    def paint_scroll_horizontal(self):
        """ Lines are fixed, columns are changing. """

        # Get visible cols
        first_visible_column, last_visible_column = self.viewer.get_visible_columns()

        if not self.viewer.color_off:
            # Right or left
            # If first_visible_column or last_visible_column has increased, it has moved to the right.
            if first_visible_column > self.prev_first_visible_column or last_visible_column > self.prev_last_visible_column: # right
                uncolored_col_range = range(max(self.prev_last_visible_column, first_visible_column), last_visible_column+1)
            elif first_visible_column < self.prev_first_visible_column or last_visible_column < self.prev_last_visible_column: # left
                uncolored_col_range = range(first_visible_column, min(self.prev_first_visible_column, last_visible_column)+1)
            else:
                return # If both are the same, there's no need to color. The shift is very slight, and no new letter has appeared.

            no_of_line = int(self.viewer.seq_text_widget.index('end-1c').split('.')[0])  # last empty included

            for line in range(1, no_of_line):

                for tag in self.viewer.all_color_tags:
                    self.viewer.seq_text_widget.tag_remove(tag, f"{line}.{uncolored_col_range[0]}",
                                                    f"{line}.{uncolored_col_range[-1] + 1}")

                # add new color tags
                seq_segment = self.viewer.seq_text_widget.get(f"{line}.{uncolored_col_range[0]}",
                                                              f"{line}.{uncolored_col_range[-1] + 1}")
                self.viewer.add_color_tag(seq_segment, line, uncolored_col_range[0],
                                          uncolored_col_range[-1])

            # raise TEXT_SEL to top
            self.viewer.seq_text_widget.tag_raise(self.viewer.TEXT_SEL)

        # update prevs
        self.prev_first_visible_column = first_visible_column
        self.prev_last_visible_column = last_visible_column

    def paint_resize_horizontal(self):
        """ Runs only if screen enlarges."""

        # Get visible cols
        first_visible_column, last_visible_column = self.viewer.get_visible_columns()

        if not self.viewer.color_off:
            # Right or left
            if first_visible_column == self.prev_first_visible_column and last_visible_column == self.prev_last_visible_column:
                return  # If both are the same, there's no need to paint. The shift is very slight, and no new letter has appeared.
            elif first_visible_column == self.prev_first_visible_column:  # rigth
                uncolored_column_range = range(self.prev_last_visible_column, last_visible_column + 1)
            elif last_visible_column == self.prev_last_visible_column:  # left
                uncolored_column_range = range(first_visible_column, self.prev_first_visible_column + 1)
            else:
                uncolored_column_range = None  # It won't enter here; this function will trigger when expanding from one side.

            no_of_line = int(self.viewer.seq_text_widget.index('end-1c').split('.')[0])  # last empty included

            for line in range(no_of_line):

                for tag in self.viewer.all_color_tags:
                    self.viewer.seq_text_widget.tag_remove(tag, f"{line}.{uncolored_column_range[0]}",
                                                    f"{line}.{uncolored_column_range[-1] + 1}")

                # add new color tags
                seq_segment = self.viewer.seq_text_widget.get(f"{line}.{uncolored_column_range[0]}",
                                                              f"{line}.{uncolored_column_range[-1] + 1}")
                self.viewer.add_color_tag(seq_segment, line, uncolored_column_range[0],
                                          uncolored_column_range[-1])

            # raise TEXT_SEL to top
            self.viewer.seq_text_widget.tag_raise(self.viewer.TEXT_SEL)

        # update prevs
        self.prev_first_visible_column = first_visible_column
        self.prev_last_visible_column = last_visible_column

    def paint_visible_lines_all(self):
        """ Finds visible lines and columns, paint only them.  """

        # Get visible lines and rows
        _, _, first_visible_column, last_visible_column = self.viewer.get_visible_region()
        # print(first_visible_line, last_visible_line, first_visible_column, last_visible_column)
        if not self.viewer.color_off:
            total_lines = int(self.viewer.seq_names_text_widget.index('end-1c').split('.')[0])

            for line in range(1, total_lines):  # paint all lines
                # delete old color tags
                for tag in self.viewer.all_color_tags:
                    self.viewer.seq_text_widget.tag_remove(tag, f"{line}.{first_visible_column}",
                                                    f"{line}.{last_visible_column + 1}")

                # add new color tags
                seq_segment = self.viewer.seq_text_widget.get(f"{line}.{first_visible_column}",
                                                       f"{line}.{last_visible_column + 1}")
                self.viewer.add_color_tag(seq_segment, line, first_visible_column, last_visible_column)

                # self.viewer.seq_text_widget.update_idletasks()  # Doing it this way paints row by row.

            # raise TEXT_SEL to top (selected be rendered above color)
            self.viewer.seq_text_widget.tag_raise(self.viewer.TEXT_SEL)

        # update prevs
        self.prev_first_visible_column = first_visible_column
        self.prev_last_visible_column = last_visible_column

    def yview_text(self, *args):
        end_time = time.time()
        runtime = end_time - self.viewer.start_time
        print("runtime", runtime)
        self.viewer.start_time = time.time()

        # Note: vertical scroll is simulated in moveto and moveto_small methods. Visible
        # print(args)
        if args[0] == 'moveto':
            f = float(args[1])
            if f == self.last_f:  # skip if same f
                # print("same f")
                return

            self.moveto_pos(f)

        elif args[0] == 'scroll':
            n = int(args[1])
            what = args[2]
            if what == tk.PAGES:
                self.scroll_pages(n)
            elif what == tk.UNITS:
                self.scroll_units(n)

    def set(self, a, b):
        """
        # Normally, the text widget sets the scrollbar based on the new position, but in our case, the tw doesn't contain all the sequences.
        # That's why we calculate and set the scrollbar ourselves, in the set_scrollbar function.

        """
        # print(a, b)
        # self.viewer.v_scrollbar_seq_names.set(a,b)
        pass

    def moveto_pos(self, f):
        """
        Move to pos. in 2 ways:
        - delete all and insert from scratch (if no intersection with previous region)
        - keep intersect and add/delete seqs
        :param f:
        :return:
        """
        # check limits
        f = self.check_f_limits(f)

        # round f
        # print("before",f)
        # f = self.round_f(f)
        # print("after",f)

        # skip if same f
        if f == self.last_f:
            # print("same f")
            return

        # find direction
        direction = 0
        if f > self.last_f:  # aşağı gitmiş
            direction = -1
        elif f < self.last_f:  # yukarı gitmiş
            direction = +1

        # check if there is intersection
        intersect, add_range, delete = self.check_intersection(f, direction)
        # print("inter:", direction, intersect, add_range, delete)

        if intersect == 0:
            # just moveto that pos. if no intersection exist
            self.moveto(f)
        else:  # keep intersection, add and delete other lines
            self.moveto_small(f, direction, add_range, delete)

    def scroll_pages(self, n):
        # find f
        f = self.last_f + n * self.window_size

        # check limits
        f = self.check_f_limits(f)

        # just moveto that pos.
        self.moveto(f)

    def scroll_units(self, n):
        # find f
        single_line_size = 1 / self.total_line_size
        f = self.last_f + n * single_line_size

        # moveto_pos
        self.moveto_pos(f)

    def check_f_limits(self, f):
        f = max(0, f)  # to fix cases where it returns -0.00... or when calling -n page.
        f = min(f, 1-self.window_size)
        return f

    def find_first_last_lines(self, f):
        # find first line (flaat)
        first_line = f * self.total_line_size + 1
        # print("first_line", first_line)

        # find first and last lines (int)
        first_line_int = math.floor(first_line)
        offset = first_line - first_line_int
        last_line_int = min(first_line_int + self.visible_count_int - 1, self.total_seq_size)
        # print("first-last", first_line_int, last_line_int)
        return first_line_int, last_line_int, offset

    def check_intersection(self, f, direction):
        """
        :param f: new position that will be checked against last one
        :return:
        """
        # find first and last lines for new position
        first_line_int, last_line_int, offset = self.find_first_last_lines(f)

        # find intersection, lines to add, lines to remove according to last position (self. variables)
        if direction == -1:  # It has moved down.
            intersect = max(0, self.last_line_int - first_line_int + 1)
            add_range = range(self.last_line_int+1, last_line_int+1)
            delete = first_line_int - self.first_line_int
            return intersect, add_range, delete

        elif direction == +1:  # It has moved up.
            intersect = max(0, last_line_int - self.first_line_int + 1)
            add_range = range(first_line_int, self.first_line_int)
            delete = self.last_line_int - last_line_int
            return intersect, add_range, delete

        return 0, 0, 0

    def update_hit(self, f):
        # check hit bottom
        if f + self.window_size >= 1:
            self.hit_bottom = True
        else:
            self.hit_bottom = False

        # check hit top
        if f <= 0:
            self.hit_top = True
        else:
            self.hit_top = False

    def check_hit(self, f):  # check if hit bottom or top
        if (f >= self.last_f and self.hit_bottom) or (f <= self.last_f and self.hit_top):
            # print("hit end")
            return True
        return False

    def set_scrollbar(self, f):
        self.viewer.v_scrollbar_seq_names.set(f, f + self.window_size)
        # print(f, f + self.window_size)

    def adjust_offset(self, offset):
        # print("offset", offset)
        tw_pos = offset / self.viewer.get_no_of_lines(
            self.viewer.seq_text_widget)  # When reaching the bottom, we print fewer sequences than visible_count_int.
        self.viewer.seq_names_text_widget.yview_moveto(tw_pos)
        self.viewer.seq_text_widget.yview_moveto(tw_pos)

    def moveto(self, f):  # This function assumes that a valid f value is given. Always invoke "check_f_limits" first !
        # print("moveto fun")
        if self.check_hit(f):
            return

        # Clear the text widgets
        self.viewer.seq_names_text_widget.configure(state="normal")  # change to normal in order to edit text widget.
        self.viewer.seq_names_text_widget.delete("1.0", tk.END)
        self.viewer.seq_text_widget.configure(state="normal")  # change to normal in order to edit text widget.
        self.viewer.seq_text_widget.delete("1.0", tk.END)

        self.first_line_int, self.last_line_int, offset = self.find_first_last_lines(f)
        # print("lines:", self.first_line_int, self.last_line_int)

        # add seqs
        self.curr_records.clear()
        for seq_id in range(self.first_line_int, self.last_line_int + 1):
            # print("inserted seq_id", seq_id)
            seq_record = self.data.get_record(seq_id)
            self.curr_records.append(seq_record)
            self.data.add_to_cache(seq_id, seq_record)

            self.viewer.display_single_seq_name(seq_id, seq_record.id, tk.END)
            self.viewer.display_single_seq_content(seq_record.seq, tk.END)
            if not self.viewer.color_off:
                last_line = int(self.viewer.seq_text_widget.index('end-1c').split('.')[0]) - 1
                seq_segment = self.viewer.seq_text_widget.get(f"{last_line}.{self.prev_first_visible_column}",
                                                              f"{last_line}.{self.prev_last_visible_column + 1}")
                self.viewer.add_color_tag(seq_segment, last_line, self.prev_first_visible_column, self.prev_last_visible_column)

        # Update the scroll region after inserting text
        self.viewer.seq_names_text_widget.update_idletasks()
        self.viewer.seq_text_widget.update_idletasks()

        # Disable editing.
        self.viewer.seq_names_text_widget.configure(state="disabled")
        self.viewer.seq_text_widget.configure(state="disabled")

        # move tw (How much the first line will be visible)
        if f + self.window_size >= 1:  # hit bottom
            # print("adjust_offset")
            self.adjust_offset(offset)

        # set scrollbar
        self.set_scrollbar(f)

        # check if hit bottom or top
        self.update_hit(f)

        # update last pos
        self.last_f = f

        # index line
        self.viewer.display_index_line()

        # column jump
        self.check_column_jump()

    def moveto_small(self, f, direction, add_range, delete):
        # print("moveto_small fun", self.prev_first_visible_column, self.prev_last_visible_column)
        if self.check_hit(f):
            print("check_hit", f, self.last_f)
            return

        self.first_line_int, self.last_line_int, offset = self.find_first_last_lines(f)
        # print("lines:", self.first_line_int, self.last_line_int)

        # if there are lines to add or delete
        if len(add_range) > 0 or delete > 0:

            # Enable edit
            self.viewer.seq_names_text_widget.configure(state="normal")  # change to normal in order to edit text widget.
            self.viewer.seq_text_widget.configure(state="normal")  # change to normal in order to edit text widget.

            if direction == -1: # Move down
                for seq_id in add_range:  # add to bottom
                    # print("inserted seq_id", seq_id)
                    seq_record = self.data.get_record(seq_id)
                    self.curr_records.append(seq_record)
                    self.data.add_to_cache(seq_id, seq_record)

                    self.viewer.display_single_seq_name(seq_id, seq_record.id, tk.END)
                    self.viewer.display_single_seq_content(seq_record.seq, tk.END)

                    if not self.viewer.color_off:
                        last_line = int(self.viewer.seq_text_widget.index('end-1c').split('.')[0]) - 1

                        # Get visible cols
                        seq_segment = self.viewer.seq_text_widget.get(f"{last_line}.{self.prev_first_visible_column}",
                                                                      f"{last_line}.{self.prev_last_visible_column + 1}")
                        self.viewer.add_color_tag(seq_segment, last_line, self.prev_first_visible_column,
                                                  self.prev_last_visible_column)

                for i in range(delete):  # delete from top
                    # print("deleted seq_id", i)
                    self.curr_records.pop(0)

                    self.viewer.delete_line(self.viewer.seq_names_text_widget, 1)
                    self.viewer.delete_line(self.viewer.seq_text_widget, 1)

            elif direction == 1: # Move up
                last_line = self.viewer.get_no_of_lines(self.viewer.seq_text_widget)-1
                for i in range(delete):  # delete from bottom
                    # print("deleted seq_id", i)
                    self.curr_records.pop(-1)

                    self.viewer.delete_line(self.viewer.seq_names_text_widget, last_line)
                    self.viewer.delete_line(self.viewer.seq_text_widget, last_line)
                    last_line -= 1

                for seq_id in reversed(add_range):  # add to top
                    # print("inserted seq_id", seq_id)
                    seq_record = self.data.get_record(seq_id)
                    self.curr_records.insert(0, seq_record)
                    self.data.add_to_cache(seq_id, seq_record)

                    self.viewer.display_single_seq_name(seq_id, seq_record.id, "1.0")
                    self.viewer.display_single_seq_content(seq_record.seq, "1.0")
                    if not self.viewer.color_off:
                        # Get visible cols
                        seq_segment = self.viewer.seq_text_widget.get(f"{1}.{self.prev_first_visible_column}",
                                                                      f"{1}.{self.prev_last_visible_column + 1}")
                        self.viewer.add_color_tag(seq_segment, 1, self.prev_first_visible_column, self.prev_last_visible_column)

            # Update the scroll region after inserting text
            self.viewer.seq_names_text_widget.update_idletasks()
            self.viewer.seq_text_widget.update_idletasks()

            # Disable editing.
            self.viewer.seq_names_text_widget.configure(state="disabled")
            self.viewer.seq_text_widget.configure(state="disabled")

        # check hit bottom
        if f + self.window_size >= 1:
            # print("adjust_offset")
            self.adjust_offset(offset)  # Move the text widget (how much the first line will be visible)

        # set scrollbar
        self.set_scrollbar(f)

        # check if hit bottom or top
        self.update_hit(f)
        # print("hit_bottom", self.hit_bottom)

        # update last pos
        self.last_f = f

        # index line
        self.viewer.display_index_line()

        # column jump
        self.check_column_jump()

    def check_column_jump(self):
        """
        Normally visible columns do not change while scrolling vertically. But if the last visible column does not exist
        in next visible area (occurs when last visible column is beyond the new sequences length ), the visible column
        range jump to left.  This cause some already added seqs to have uncolored letters. We must paint them.
        :return:
        """
        last_visible_column = int(
            self.viewer.seq_ruler_widget.index(f"@0,{self.viewer.seq_ruler_widget.winfo_width()}").split('.')[1])
        if last_visible_column != self.prev_last_visible_column:
            print("check_column_jump")
            # self.prev_last_visible_column = last_visible_column
            # self.prev_first_visible_column = int(self.viewer.seq_ruler_widget.index("@0,0").split('.')[1])
            self.paint_visible_lines_all()

    def add_colored_seq(self, seq, pos):
        # Loop through the sequence and group contiguous identical letters, then add corr tag
        last_letter = None
        group_start = 0
        inserted_len = 0
        for col_index, letter in enumerate(seq):
            if letter != last_letter:
                if last_letter is not None:
                    # Apply the tag for the previous group of identical letters
                    tag_name = f"tag_{last_letter.upper()}"
                    self.viewer.seq_text_widget.insert(f"{pos}+{inserted_len}c", seq[group_start:col_index], tag_name)
                    inserted_len += col_index-group_start
                last_letter = letter
                group_start = col_index

        # Apply the tag for the last group of identical letters
        if last_letter is not None:
            tag_name = f"tag_{last_letter.upper()}"
            self.viewer.seq_text_widget.insert(f"{pos}+{inserted_len}c", seq[group_start:], tag_name)
            inserted_len += len(seq) - group_start

        # add new line char
        self.viewer.seq_text_widget.insert(f"{pos}+{inserted_len}c", "\n", ())  # insert without tag

    def height_decreased(self):
        print("height_decreased")
        self.update_state()
        direction = +1  # If it has moved up, we will delete from the bottom.
        self.hit_bottom = False

        # check if there is intersection
        intersect, add_range, delete = self.check_intersection(self.last_f, direction)
        # print("inter:", direction, intersect, add_range, delete)

        # keep intersection, add and delete other lines
        self.moveto_small(self.last_f, direction, add_range, delete)

    def height_increased(self):
        # print("height_increased")
        self.update_state()

        if self.hit_bottom:
            self.hit_bottom = False
            direction = +1  # yukarı çıkmış gibi yukarı eklicez
            f = 1-self.window_size
        else:
            self.hit_top = False
            direction = -1  # aşağı inmiş gibi aşağı eklicez
            f = self.last_f

        # print("f", f)

        # check if there is intersection
        intersect, add_range, delete = self.check_intersection(f, direction)
        # print("inter:", direction, intersect, add_range, delete)

        # keep intersection, add and delete other lines
        self.moveto_small(f, direction, add_range, delete)

    def width_increased(self):
        # print("width_increased")
        self.paint_resize_horizontal()

    def width_decreased(self):
        # print("width_decreased")
        self.update_visible_columns()  # just update visible colns

    def update_visible_columns(self):
        """ Just updates visible columns.  """
        # Get visible cols
        first_visible_column, last_visible_column = self.viewer.get_visible_columns()

        # update prevs
        self.prev_first_visible_column = first_visible_column
        self.prev_last_visible_column = last_visible_column

    def on_resize(self, event):
        # print(f"Window resized to {event.width}x{event.height}")
        if not self.viewer.seqs_displayed:
            return

        # if event.height != self.last_tw_height or event.width != self.last_tw_width:
            # print("height", event.height, self.last_tw_height, "width", event.width, self.last_tw_width)

        # ---------SEQ ADD-DELETE-------
        # Display index line is called in moveto_small.
        # height increased
        if event.height > self.last_tw_height:
            self.height_increased()

        # height decreased
        elif event.height < self.last_tw_height:
            self.height_decreased()

        # -------- ONLY PAINT----------
        # width increased
        if event.width > self.last_tw_width:
            self.width_increased()

        # width decreased
        elif event.width < self.last_tw_width:
            self.width_decreased()

        # --------------------------
        self.last_tw_height = event.height
        self.last_tw_width = event.width

    def font_size_increased(self):  # DISABLED
        # same behaviour with height and width decrease
        # self.height_decreased() -- this is already called in on_resize after font size increased
        # self.paint_visible_lines_all()  # for columns
        pass

    def font_size_decreased(self):  # DISABLED
        # same behaviour with height and width increase
        # self.height_increased()  -- this is already called in on_resize after font size decreased
        # self.paint_visible_lines_all()  # for columns
        pass

    def set_mouse_wheel(self, text_widget):
        # Bind the scroll events
        text_widget.bind("<MouseWheel>", self.sync_scroll)

        # Bind for Linux systems (Button-4 and Button-5 for mouse wheel scroll)
        text_widget.bind("<Button-4>", self.sync_scroll_linux)
        text_widget.bind("<Button-5>", self.sync_scroll_linux)

    def sync_scroll(self, event):
        # Determine platform
        sys_platform = platform.system()

        if sys_platform == 'Windows':
            # On Windows, delta is a multiple of 120
            delta = int(-1 * (event.delta / 120))
        elif sys_platform == 'Darwin':
            # On macOS, delta is a multiple of 1
            delta = int(-1 * event.delta)
        else:
            # On Linux, use <Button-4> and <Button-5>
            return

        self.scroll_units(delta)
        return "break"

    def sync_scroll_linux(self, event):
        # On Linux, Button-4 is scroll up and Button-5 is scroll down
        if event.num == 4:
            self.scroll_units(-1)
        elif event.num == 5:
            self.scroll_units(+1)
        return "break"

    def update_color_mode(self):
        self.hit_top = False
        self.hit_bottom = False
        self.moveto(self.last_f)  # print again

    def show_mouse_position(self, event):
        # Get the line and column from the mouse position
        line_col = self.viewer.seq_text_widget.index(f"@{event.x},{event.y}")
        line = int(line_col.split(".")[0]) + self.first_line_int - 1
        col = line_col.split(".")[1]
        # print(int(line_col.split(".")[0]), self.first_line_int)

        # Print the line and column in the format "line.column"
        # print(f"Mouse is over line.column: {line_col}")
        self.viewer.coordinate_var.set(f"Row: {line} Column: {int(col) + 1}  ")

    def navigate_to_pos(self, line, column):
        f = (line - 1)/self.total_line_size  # find f when line is at the top of window
        f = f - self.window_size/2  # find f when line is at middle point of screen

        # check limits and go
        f = self.check_f_limits(f)
        self.moveto(f)

        rel_line = line - (self.first_line_int - 1)  # line number relative to tw
        position = self.viewer.seq_text_widget.index(f"{rel_line}.{column}")  # bound dışındakiler bound olarak dönüyor
        self.viewer.see_position(*self.viewer.index_to_line_col(position))

        self.paint_visible_lines_all()

    def replace_selected_text(self, replacement_char, check_editor_mode=True):
        if check_editor_mode and not self.viewer.editor_mode_var.get():
            print("Editor mode off.")
            return

        if not self.viewer.any_selection(self.viewer.seq_text_widget):
            print("No selected text.")
            return

        # Get all selected ranges in text widget (for re-select replaced region)
        selection = self.viewer.get_selection(self.viewer.seq_text_widget)

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.viewer.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        self.viewer.seq_text_widget.configure(state="normal")

        # Iterate over each selected line and delete
        for start, end in self.viewer.get_selection_generator(self.viewer.seq_text_widget, tag_ranges=selection):
            # count no of letters to replace (same for each line)
            count = len(self.viewer.seq_text_widget.get(start, end))

            # Delete the characters within the tagged range
            # self.seq_text_widget.tag_remove("color", start, end)
            self.viewer.seq_text_widget.delete(start, end)

            # Insert the replacement character at the start index COLORED INDEX
            tag_name = f"tag_{replacement_char.upper()}"
            self.viewer.seq_text_widget.insert(start, count * replacement_char, tag_name)

            # add to edited dict
            edited_line, _ = self.viewer.index_to_line_col(start)
            self.add_to_edited(edited_line)

        self.viewer.seq_text_widget.configure(state="disabled")

        # reselect replaced text
        for start, end in self.viewer.get_selection_generator(self.viewer.seq_text_widget, tag_ranges=selection):
            self.viewer.seq_text_widget.tag_add(self.viewer.TEXT_SEL, start, end)

        # self.seq_text_widget.update_idletasks()
        self.viewer.master.unsaved_file = True

    def add_to_edited(self, edited_line):
        seq_id = self.viewer.get_seq_id(edited_line)
        seq_record = self.curr_records[edited_line - 1]
        seq_record.letter_annotations = {}
        seq_record.seq = self.viewer.get_line(self.viewer.seq_text_widget, edited_line)
        self.data.add_to_edited(seq_id, seq_record)

    def delete_selected_letters(self):
        # Get all selected ranges in text widget
        selection = self.viewer.get_selection(self.viewer.seq_text_widget)
        if not selection:
            print("No selected text.")
            return

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.viewer.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        # Get all selected ranges in text widget
        fully_deleted_lines = self.viewer.check_whole_line_deleted(selection)

        if len(fully_deleted_lines) == 0:
            self.viewer.seq_text_widget.configure(state="normal")

            # Iterate over each tagged range
            for start, end in self.viewer.get_selection_generator(self.viewer.seq_text_widget, selection):
                # Delete the characters within the tagged range
                self.viewer.seq_text_widget.delete(start, end)

                # add to edited dict
                edited_line, _ = self.viewer.index_to_line_col(start)
                self.add_to_edited(edited_line)

            self.viewer.seq_text_widget.configure(state="disabled")

            self.viewer.master.unsaved_file = True
            self.viewer.display_index_line()
            self.paint_visible_lines_all()

            self.viewer.remove_selection(self.viewer.seq_names_text_widget) # seq_text already deleted, ruler cannot be selected
        else:
            print("Cannot delete whole sequence.")
            return

    def insert_unknown_char(self):
        """
          We insert only if there is selected text. So, seq_text_widget will always be edited. seq_names_text_widget not.
         :return:
         """
        if not self.viewer.any_selection(self.viewer.seq_text_widget):
            return

        selection = self.viewer.get_selection(self.viewer.seq_text_widget)

        # Check if letter annot exist in selected seqs, if yes ask user to delete annot
        confirm_edit = self.viewer.check_seletion_annot(selection)
        if not confirm_edit:
            print("User denied deleting letter annots. Replace was not applied")
            return

        self.viewer.seq_text_widget.configure(state="normal")

        # add indel (?) at start of selected block
        for start, end in self.viewer.get_selection_generator(self.viewer.seq_text_widget, selection):
            line_number, _ = self.viewer.index_to_line_col(start)
            self.viewer.seq_text_widget.insert(start, "?")

            # add to edited dict
            self.add_to_edited(line_number)

        self.viewer.seq_text_widget.configure(state="disabled")

        self.viewer.remove_selection(self.viewer.seq_text_widget)
        self.viewer.remove_selection(self.viewer.seq_ruler_widget)
        self.viewer.remove_selection(self.viewer.seq_names_text_widget)

        self.viewer.master.unsaved_file = True
        self.viewer.display_index_line()
        self.paint_visible_lines_all()

    def open_annot(self):
        if not self.viewer.any_selection(self.viewer.seq_names_text_widget):
            print("no selected viewer.seq for annot")
            return

        # get first selected seq
        selection = self.viewer.get_selection(self.viewer.seq_names_text_widget)
        start = selection[0]  # start of first line
        line_number, _ = self.viewer.index_to_line_col(start)
        seq_id = self.viewer.get_seq_id(line_number)

        # Get seq_record
        print("Annot Seq ID :", seq_id)
        seq_record = self.data.get_record(seq_id)

        # check if seq is different in text widget
        line_text = self.viewer.get_line(self.viewer.seq_text_widget, line_number)
        show_letter_annot = (line_text == str(seq_record.seq))
        print("show_letter_annot:", seq_id, show_letter_annot)

        # we give the "seq_record" object to annot window and it updates seq_record if user click "save and close".
        AnnotationViewer(seq_id, seq_record, line_text, self.viewer.master, show_letter_annot, self.annot_seq_id_edit,
                         self.save_current_workspace, line_number)

    def annot_seq_id_edit(self, line_number, seq_id, seq_record_id, seq_record_seq):
        """
        If seq id is edited in annot window, left panel is updated.
        :param line_number:
        :param seq_id:
        :param seq_record_id:
        :return:
        """
        self.viewer.seq_names_text_widget.configure(state="normal")

        # Delete old line
        self.viewer.delete_line(self.viewer.seq_names_text_widget, line_number)
        self.viewer.delete_line(self.viewer.seq_text_widget, line_number)

        # Insert new line
        position = f"{line_number}.0"
        self.viewer.display_single_seq_name(seq_id, seq_record_id, position)
        self.viewer.display_single_seq_content(seq_record_seq, position)

        self.viewer.seq_names_text_widget.configure(state="disabled")
        self.viewer.master.unsaved_file = True

    def save_current_workspace(self):
        # undo redo does not exist in index mode, so no need to save workspace
        pass

    def move_selection(self, direction):
        if not self.viewer.any_selection(self.viewer.seq_text_widget):
            return

        _, _ ,_ , start_index = self.viewer.get_selection_bounds()
        line_num, col_num = map(int, start_index.split("."))
        col_num -= 1

        # Move the tag position depending on the arrow key
        if direction == "Left":
            col_num -= 1
        elif direction == "Right":
            col_num += 1
        elif direction == "Up":
            line_num -= 1
        elif direction == "Down":
            line_num += 1

        # Construct the new position index
        new_pos = f"{line_num}.{col_num}"
        # print(new_pos)

        # Check if the new position is within valid bounds
        letter = self.viewer.seq_text_widget.get(new_pos)
        if line_num > 0 and col_num >= 0 and letter != "\n" and letter != "":

            if direction in ["Up", "Down"]:
                last_full_visible_line = math.floor(self.visible_count)
                if direction == "Down" and line_num == last_full_visible_line + 1:
                    return  # no move to outside of visible area
                else:
                    self.viewer.highlight_pos(line_num, col_num, highlight_len=1)

            else:
                self.viewer.see_position(line_num, col_num, highlight_len=1)
                self.paint_scroll_horizontal()


class TranslateOptions(ModalWindow):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(self.parent, resizable=False)
        self.title("Translate Options")

        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(3, weight=1)
        self.grid_columnconfigure(0, weight=1, uniform='a')
        self.grid_columnconfigure(1, weight=1, uniform='a')
        self.grid_columnconfigure(2, weight=1, uniform='a')
        # self.grid_columnconfigure(3, weight=1, uniform='a')
        # self.grid_columnconfigure(4, weight=1, uniform='a')

        self.genetic_codes = ['Standard', 'SGC0', 'Vertebrate Mitochondrial', 'SGC1', 'Yeast Mitochondrial', 'SGC2', 'Mold Mitochondrial',
         'Protozoan Mitochondrial', 'Coelenterate Mitochondrial', 'Mycoplasma', 'Spiroplasma', 'SGC3',
         'Invertebrate Mitochondrial', 'SGC4', 'Ciliate Nuclear', 'Dasycladacean Nuclear', 'Hexamita Nuclear', 'SGC5',
         'Echinoderm Mitochondrial', 'Flatworm Mitochondrial', 'SGC8', 'Euplotid Nuclear', 'SGC9', 'Bacterial',
         'Archaeal', 'Plant Plastid', 'Alternative Yeast Nuclear', 'Ascidian Mitochondrial',
         'Alternative Flatworm Mitochondrial', 'Blepharisma Macronuclear', 'Chlorophycean Mitochondrial',
         'Trematode Mitochondrial', 'Scenedesmus obliquus Mitochondrial', 'Thraustochytrium Mitochondrial',
         'Pterobranchia Mitochondrial', 'Candidate Division SR1', 'Gracilibacteria', 'Pachysolen tannophilus Nuclear',
         'Karyorelict Nuclear', 'Condylostoma Nuclear', 'Mesodinium Nuclear', 'Peritrich Nuclear',
         'Blastocrithidia Nuclear', 'Balanophoraceae Plastid', 'Cephalodiscidae Mitochondrial']

        # init
        self.is_closed_succesfully = False
        self.genetic_code = self.genetic_codes[0]
        self.to_stop = False
        self.cds = False

        # Create and place labels and entry widgets
        tk.Label(self, text="Genetic Code:").grid(row=0, column=0, padx=5, pady=10, sticky=tk.W)

        self.genetic_code_var = tk.StringVar(self)
        self.genetic_code_var.set(self.genetic_code)
        self.genetic_code_menu = tk.OptionMenu(self, self.genetic_code_var, *self.genetic_codes)
        self.genetic_code_menu.grid(row=0, column=1, columnspan=2, padx=5, pady=10, sticky=tk.W)

        # Checkboxes
        self.to_stop_var = tk.BooleanVar()
        self.cds_var = tk.BooleanVar()

        self.to_stop_check = ttk.Checkbutton(self, text="Terminate at the first in frame stop codon", variable=self.to_stop_var)
        self.cds_check = ttk.Checkbutton(self, text="Complete Coding Sequence (CDS)", variable=self.cds_var)

        self.to_stop_check.grid(row=1, column=0, columnspan=3, sticky="w", padx=5, pady=2)
        self.cds_check.grid(row=2, column=0, columnspan=3, sticky="w", padx=5, pady=2)

        # Create and place the convert button
        tk.Button(self, text="Confirm", command=self.confirm).grid(row=3, column=1, padx=5, pady=20, sticky="n")

        # Pause execution in the main program until this window is closed
        self.parent.wait_window(self)

    def confirm(self):
        self.genetic_code = self.genetic_code_var.get()
        self.to_stop = self.to_stop_var.get()
        self.cds = self.cds_var.get()
        self.is_closed_succesfully = True
        self.destroy()
