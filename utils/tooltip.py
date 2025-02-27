import tkinter as tk


class Tooltip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tooltip_window = None
        self.active = True
        widget.bind("<Enter>", self.show_tooltip)
        widget.bind("<Leave>", self.hide_tooltip)

    def show_tooltip(self, event):
        if self.active:
            x = event.x_root + 20
            y = event.y_root + 20
            self.tooltip_window = tw = tk.Toplevel(self.widget)
            tw.transient(self.widget)
            tw.wm_overrideredirect(True)
            tw.wm_geometry(f"+{x}+{y}")
            label = tk.Label(tw, text=self.text, background="#ffffe0", relief='solid', borderwidth=1)
            label.pack()

    def hide_tooltip(self, event):
        if self.tooltip_window:
            self.tooltip_window.destroy()
            self.tooltip_window = None

    def activate(self):
        self.active = True

    def deactivate(self):
        self.active = False
        self.hide_tooltip(None)  # Ensure the tooltip is hidden when deactivated

    def change_text(self, text):
        self.text = text
        if self.tooltip_window:
            label = self.tooltip_window.winfo_children()[0]
            label.config(text=text)
