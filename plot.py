import time, sys, os, Tkinter, numpy, math
from subprocess import Popen, PIPE

class Plot:
    """
    Plot a vector or list of vectors. Assumes that all vectors in the list
    are the same type (Float or Complex) Prompts for which ones to show.
    """
    def __init__(self, data, **kwargs):
        self.quiet = kwargs.get('quiet', True)
        self.commands = kwargs.get('commands', '')
        self.linewidth=kwargs.get('linewidth', 1.0)
        self.style = kwargs.get('style', '')
        self.args = kwargs
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
        if len(self.data) > 0:
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
        self.gnuplot = Popen([
            'export DYLD_LIBRARY_PATH= ;export LD_LIBRARY_PATH= ;gnuplot',
            '-geometry 1200x1000+200+0'],
                             shell=True, stdin=PIPE)
    
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
        self.figure = MF = MatplotFigure(add_subplot=False)
        MF.axis = MF.figure.add_axes( [0.07, 0.07, 0.8, 0.9] )
        n = len(self.data)
        self.funcs_to_show = [Tk.BooleanVar(value=True) for i in range(n)]
        func_selector_frame = ttk.Frame(MF.window)
        for i in range(n):
            var = self.funcs_to_show[i]
            button = ttk.Checkbutton(func_selector_frame, text='%d'% i,
                                     variable=var, command=self.show_plots)
            button.grid(column=0, row=i, sticky=(Tk.N, Tk.W))
            #frame = ttk.Frame(func_selector_frame, 
        func_selector_frame.grid(column=1, row=0, sticky=(Tk.N))
        MF.window.columnconfigure(1, weight=0)

    def test(self):
        return [v.get() for v in self.funcs_to_show]

    def color(self, i):
        from matplotlib.cm import gnuplot2
        return gnuplot2(i/8.0 - math.floor(i/8.0))

    def split_data(self, data):
        """
        The data has None entries between points which should not
        be connected by arcs in the picture.  This method splits the
        data at the None entries, and builds the x and y lists for the
        plotter.
        """
        result = []
        if self.type == 'complex':
            x_list, y_list = [], [] 
            for d in data:
                if d is None and len(x_list) > 1:
                    result.append( (x_list, y_list) )
                    x_list, y_list = [], [] 
                else:
                    x_list.append(d.real)
                    y_list.append(d.imag)
        else:
            x_list, y_list = [], [] 
            for n, d in enumerate(data):
                if d is None and len(x_list) > 1:
                    result.append( (x_list, y_list) )
                    x_list, y_list =[], []
                else:
                    x_list.append(n)
                    y_list.append(d)
        result.append( (x_list, y_list) )
        return result
                    
    def create_plot(self):
        axis = self.figure.axis
        axis.clear()
        for i, data in enumerate(self.data):
            if self.funcs_to_show[i].get():
                lists = self.split_data(data)
                for X, Y in lists:
                    axis.plot(X, Y, color=self.color(i),
                              linewidth=self.linewidth, label='%d' % i)

        # Configure the plot based on keyword arguments
        limits = self.args.get('limits', None)
        xlim, ylim = limits if limits else (axis.get_xlim(), axis.get_ylim())

        margin_x, margin_y = self.args.get('margins', (0.1, 0.1))
        sx = ( xlim[1] - xlim[0])*margin_x
        xlim = (xlim[0] - sx, xlim[1] + sx)
        sy = (ylim[1] - ylim[0])*margin_y
        ylim = (ylim[0] - sy, ylim[1] + sy)
        axis.set_xlim(*xlim)
        axis.set_ylim(*ylim)

        axis.set_aspect(self.args.get('aspect', 'auto'))
        legend = axis.legend(loc='upper left', bbox_to_anchor = (1.0, 1.0))
        decorator = self.args.get('decorator', None)

        title = self.args.get('title', None)
        if title:
            self.figure.window.title(title)

        self.figure.draw()

    def show_plots(self):
        self.create_plot()
    

if __name__ == "__main__":
    MyPlot = MatplotPlot
    #float_data = numpy.random.random( (10,) )
    #MyPlot(float_data)
    P = MyPlot([[ a+b*1j for a, b in numpy.random.random( (10,2) )] for i in range(15)])
