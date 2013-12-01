"""
Plotting using matplotlib (included in Sage) and Tkinter.
"""

# Load Tkinter
import sys, os
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk
import ttk

# Load MatplotLib
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg


class MatplotFigure:
    def __init__(self, add_subplot=True, root=None):
        figure = matplotlib.figure.Figure(figsize=(6,6), dpi=100)
        figure.set_facecolor('white')
        axis = figure.add_subplot(111) if add_subplot else None
        self.figure, self.axis = figure, axis
        
        window = Tk.Tk() if root is None else Tk.Toplevel(root)
        figure_frame = ttk.Frame(window)
        canvas = FigureCanvasTkAgg(figure, master=figure_frame)
        canvas._tkcanvas.config(highlightthickness=0, width=600, height=600)
        toolbar = NavigationToolbar2TkAgg(canvas, figure_frame)
        toolbar.pack(side=Tk.TOP, fill=Tk.X)
        canvas._tkcanvas.pack(side=Tk.TOP,  fill=Tk.BOTH, expand=1)
        toolbar.update()
        
        figure_frame.grid(column=0, row=0, sticky=(Tk.N, Tk.S, Tk.E, Tk.W))
        window.columnconfigure(0, weight=1)
        window.rowconfigure(0, weight=1)
        self.window, self.canvas, self.toolbar = window, canvas, toolbar
        self.figure_frame = figure_frame

    def draw(self):
        self.canvas.draw()

    def clear(self):
        self.axis.clear()
        self.draw()
        

if __name__ == "__main__":
    from numpy import arange, sin, pi
    MF = MatplotFigure()
    t = arange(0.0,3.0,0.01)
    s = sin(2*pi*t)
    MF.axis.plot(t,s)
    Tk.mainloop()
