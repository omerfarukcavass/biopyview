import tkinter as tk
from typing import Union


class ModalWindow(tk.Toplevel):
    """
    Use this class for all modal windows. It provides modal window specs:
    - Attach this window to its parent (transient)
    - Remove interation of user with its main window without closing (grab_set)
    - Disable/enable menu of its parent window as needed.
    - It has disable/enable menu methods  (for its children's use)
    """
    def __init__(self, parent: Union[tk.Tk, tk.Toplevel], resizable=True, min_size=None, grab_set=True):
        self.parent = parent # parent window of this modal window

        # Check if the parent has the required methods
        if not hasattr(parent, 'disable_menu') or not hasattr(parent, 'enable_menu'):
            raise TypeError("The parent must implement the 'disable_menu' and 'enable_menu'.")

        # main setup
        super().__init__(self.parent)

        # When a window is made transient to another window using the transient() method,
        # it means that it will appear "attached" to the parent window and
        # will be treated as a child window in terms of focus and stacking order.
        # However, it does not stop the execution of code in the parent window.
        # To stop the execution of code in the parent window, use:
        #       self.parent.wait_window(self)
        #
        self.transient(self.parent)

        # When a window is grabbed using grab_set(), it can have a modal behavior,
        # meaning that it prevents the user from interacting with other windows
        # until the grabbed window is closed or grab_release() is called.
        #
        # However, calling grab_set() itself doesn't directly affect the execution
        # of code in the parent window. The code in the parent window will continue
        # to execute as normal unless explicitly blocked by waiting for the grabbed
        # window to be destroyed or released. That is, the code block below creating
        # modal window will run after modal window is created.
        # To stop the execution of code in the parent window, use:
        #       self.parent.wait_window(self)
        #
        if grab_set:
            self.grab_set()
            self.is_grab_set = True
        else:
            self.is_grab_set = False

        # Disable window resizing
        if isinstance(resizable, tuple):
            self.resizable(*resizable)
        else:
            self.resizable(resizable, resizable)

        # min size
        if min_size is not None:
            self.minsize(min_size[0], min_size[1])

        self.bind("<Destroy>", self.on_destroy)

        # Disable parent window menu
        self.parent.disable_menu()  # checked above

    def on_destroy(self, event):
        # This function is called ONCE FOR EACH WIDGET if window is closed. It is called whenever a widget is destroyed,
        # in window. So we check below conditon to run only for Toplevel window not its widgets.
        if event.widget == self:

            # Enable parent window menu before closing
            self.parent.enable_menu()  # checked above

            # If parent is also a modal window and parent DID call grab_set, make parent grab the events
            if isinstance(self.parent, ModalWindow) and self.parent.is_grab_set:
                self.parent.grab_set()

