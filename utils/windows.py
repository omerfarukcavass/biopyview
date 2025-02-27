from utils.model_window import ModalWindow
import tkinter as tk
from typing import Callable
from typing import Union


class MessageBox(ModalWindow):
    def __init__(self, parent, text, window_title=None, block_parent=True, ok_button=True):
        self.parent = parent
        self.text = text
        super().__init__(self.parent, resizable=False)

        if window_title is not None:
            self.title(window_title)
        else:
            self.title("Error")

        # Disable interaction with other windows
        self.grab_set()

        # Add a label with the error message
        self.label = tk.Label(self, text= self.text)
        self.label.pack(pady=10, padx=10)

        # Add an OK button to close the popup
        if ok_button:
            button = tk.Button(self, text="OK", command=self.destroy)
            button.pack(pady=5)

        # Pause execution in the main program until this window is closed
        if block_parent:
            self.parent.wait_window(self)

        # Set the protocol to intercept the close event (MAKE USER CANNOT CLOSE WINDOW)
        self.protocol("WM_DELETE_WINDOW", self.on_closing)

    def update_text(self, text):
        self.label.config(text=text)

    def on_closing(self):
        # Do nothing when the close button is pressed
        pass


class TKMessageBox:
    def __init__(self, title, message, tk_messagebox: Callable, window: Union[tk.Tk, tk.Toplevel]):
        self.title = title
        self.message = message
        self.tk_messagebox = tk_messagebox

        # Check if the window has the required methods
        if not hasattr(window, 'disable_menu') or not hasattr(window, 'enable_menu'):
            raise TypeError("The parent must implement the 'disable_menu' and 'enable_menu'.")
        else:
            self.window = window

    def show(self):
        self.window.disable_menu()
        answer = self.tk_messagebox(self.title, self.message)
        self.window.enable_menu()
        return answer





