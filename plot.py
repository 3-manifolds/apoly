import time, sys, os, Tkinter, numpy
from subprocess import Popen, PIPE

class Plot:
    """
    Plot a vector or list of vectors. Assumes that all vectors in the list
    are the same type (Float or Complex) Prompts for which ones to show.
    """
    def __init__(self, data, quiet=True, commands='', linewidth=2, style='lines'):
        self.quiet = quiet
        self.commands = commands
        self.linewidth=linewidth
        self.style = style
        if isinstance(data[0], list) or isinstance(data[0], numpy.ndarray):
            self.data = data
        else:
            self.data = [data]
        duck = self.data[0][0]
        self.type = type(duck)
        if 'complex' in str(self.type):
            self.type = 'complex'
        elif 'float' in str(self.type):
            self.type = 'float'
        else:
            print 'Type is:', self.type

        self.start_plotter()
        if len(self.data) > 1:
            self.show_plots()
        else:
            self.create_plot([0])
            time.sleep(1)
        
    def __repr__(self):
        return ''

    def show_plots(self):
        if not self.quiet:
            print 'There are %d functions.'%len(self.data)
            print 'Which ones do you want to see?'
        else:
            self.create_plot( range(len(self.data)) )
        while 1:
            try:
                stuff = raw_input('plot> ')
                items = stuff.split()
                if len(items) and items[0] == 'all':
                    funcs = range(len(self.data))
                else:
                    funcs = [int(item)%len(self.data) for item in items]
                if len(funcs) == 0:
                    break
            except ValueError:
                break
            print funcs
            self.create_plot(funcs)
        return

    def start_plotter(self):
        """
        Stub for starting up the plotting window. 
        """

    def create_plot(self, funcs):
        """
        Stub for drawing the plot itself.
        """


class GnuplotPlot(Plot):
    """
    Uses gnuplot to plot a vector or list of vectors.
    Assumes that all vectors in the list are the same type (Float or Complex)
    Prompts for which ones to show.
    """
    def start_plotter(self):
        self.gnuplot = Popen(['export DYLD_LIBRARY_PATH= ; gnuplot',
                              '-geometry 1200x1000+200+0'],
                             shell=True,
                             stdin=PIPE)
    
    def create_plot(self, funcs):
        spec = []
        for n in funcs:
            spec.append('"-" t "%d" w %s lw %s'%(n, self.style, self.linewidth))
        gnuplot_input = self.commands + 'plot ' + ', '.join(spec) + '\n'
        if self.type == 'complex':
            for n in funcs:
                gnuplot_input += '\n'.join([
                    '%f %f'%(point.real, point.imag) if point is not None else ''
                    for point in self.data[n]] + ['e\n']) 
        elif self.type == 'float':
            for n in funcs:
                gnuplot_input += '\n'.join(
                    ['%f'%float(point) for point in self.data[n]] + ['e\n']) 
        else:
            print self.type
            print self.data[0]
            raise ValueError('Data must consist of vectors of real or complex numbers.')
        self.gnuplot.stdin.write(gnuplot_input)
        

class SagePlot(Plot):
    def create_plot(self, funcs):
        from sage.all import Graphics, line, colormaps, floor
        cm = colormaps['gnuplot2']

        G = Graphics()
        for f in funcs:
            if self.type == 'complex':
                points = [(d.real, d.imag) for d in self.data[f]]
            else:
                points = [ (i,d) for i, d in enumerate(self.data[f])]
            G += line(points, color=cm( f/8.0 - floor(f/8.0) )[:3],
                      thickness=self.linewidth, legend_label='%d' % f)
        G.show()

class MatplotPlot(Plot):
    #def __init__(self, data, **kwargs):
    #    Plot.__init__(self, data, **kwargs)
        
    def start_plotter(self):
        from tkplot import MatplotFigure, Tk, ttk
        self.figure = MF = MatplotFigure()
        n = len(self.data)
        self.funcs_to_show = [Tk.BooleanVar(value=True) for i in range(n)]
        func_selector_frame = ttk.Frame(MF.window)
        for i in range(n):
            var = self.funcs_to_show[i]
            button = ttk.Checkbutton(func_selector_frame, text='%d'% i, variable=var,
                                     command=self.show_plots)
            button.grid(column=0, row=i, sticky=(Tk.N, Tk.W))
            #frame = ttk.Frame(func_selector_frame, 


        func_selector_frame.grid(column=1, row=0, sticky=(Tk.N))
        MF.window.columnconfigure(1, weight=0)

    def test(self):
        return [v.get() for v in self.funcs_to_show]

    def color(self, i):
        from matplotlib.cm import gnuplot2
        import math
        return gnuplot2(i/8.0 - math.floor(i/8.0))
    def color(self, i):
        from matplotlib.cm import gnuplot2
        import math
        return gnuplot2(i/8.0 - math.floor(i/8.0))

    def create_plot(self):
        axis = self.figure.axis
        axis.clear()
        for i, data in enumerate(self.data):
            if self.funcs_to_show[i].get():
                if self.type == 'complex':
                    X, Y = [d.real for d in data], [d.imag for d in data]
                else:
                    X, Y = range(len(data)), data
                axis.plot(X, Y, color=self.color(i), linewidth=self.linewidth, label='%d' % i)
        axis.legend()
        self.figure.draw()

    def show_plots(self):
        self.create_plot()

        
        
            
            
    

if __name__ == "__main__":
    MyPlot = MatplotPlot
    #float_data = numpy.random.random( (10,) )
    #MyPlot(float_data)
    P = MyPlot([[ a+b*1j for a, b in numpy.random.random( (10,2) )] for i in range(15)])