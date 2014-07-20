#!/usr/bin/env python
"""
Script to pick a timewindow in SAC files and save
the start and end times as header variables
"""

import numpy as npy
import matplotlib.pyplot as plt
import argparse
from pysmo.sac.sacio	import sacfile
from matplotlib.widgets import SpanSelector, Button
from scipy.signal	import butter, lfilter

def getargs():
    hp = 0.02
    lp = 0.1
    n  = 3
    parser = argparse.ArgumentParser(description='Pick and save timewindow of one or multiple SAC files.',
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--filter', dest='filter', action="store_true",
                            help='Superimpose bandpass filtered data.')
    parser.add_argument('-H', '--hp', dest='hp', type=float, default=hp,
                            help='Highpass frequency.')
    parser.add_argument('-L', '--lp', dest='lp', type=float, default=lp,
                            help='Lowpass frequency.')
    parser.add_argument('-N', '--order', dest='n', type=int, default=n,
                            help='Filter Order.')
    parser.add_argument('-D', dest='header', type=str, nargs=2, default=['t0', 't1'],
                            help='Header names to store time window.')
    parser.add_argument(dest='sacfile', type=str, nargs='+')
    args = parser.parse_args()
    return args

def readsac(infile):
    sac = sacfile(infile, 'rw')
    y = npy.array(sac.data)
    x = npy.arange(sac.b, sac.b+sac.npts*sac.delta, sac.delta)
    return sac, x, y

def filterdata(sac, args):
    Nyq = 1 / (2 * sac.delta)
    W_hp = args.hp/Nyq
    W_lp = args.lp/Nyq
    n     = args.n
    [b,a] = butter(n, W_hp, btype='high')
    y2 = lfilter(b, a, sac.data)
    [b,a] = butter(n, W_lp, btype='low')
    y2 = lfilter(b, a, y2)
    return npy.array(y2)

def StartPlot(args):
    def plotData(sac, x, y, args):
        ax1.cla()
        ax2.cla()
        ax3.cla()
        ax1.set_title('Editing %s and %s in %s' % (t1, t2, sac.filename))
        line1, = ax1.plot(x, y, '-')
        line2, = ax2.plot(x, y, '-')
        line3, = ax3.plot(x, y, '-')
        if args.filter:
            y_filt = filterdata(sac, args)
            y_filt *= max(sac.data) / y_filt.max()
            line1_filt = ax1.plot (x, y_filt, 'r-', alpha=0.6)
            line2_filt = ax2.plot(x, y_filt, 'r-', alpha=0.6)
            line3_filt = ax3.plot(x, y_filt, 'r-', alpha=0.6)

        try:
            t1val, t2val = sac.__getattr__(t1), sac.__getattr__(t2)
            xmin, xmax = t1val, t2val
            ax1.fill([t1val,t2val,t2val,t1val], [min(y),min(y),max(y),max(y)],
            'g', alpha=0.2, edgecolor='g')
            ax2.fill([t1val,t2val,t2val,t1val], [min(y),min(y),max(y),max(y)],
            'g', alpha=0.2, edgecolor='g')
            ax3.fill([t1val,t2val,t2val,t1val], [min(y),min(y),max(y),max(y)],
            'g', alpha=0.2, edgecolor='g')
            tmpmin = t1val - 200
            tmpmax = t2val + 200
            if tmpmin < sac.b:
                tmpmin = sac.b
            if tmpmax > sac.e:
                tmpmax = sac.e
            ax2.set_xlim(tmpmin, tmpmax)
            ax3.set_xlim(t1val, t2val)
        except ValueError:
            ax2.set_xlim(sac.b, sac.e)
            ax3.set_xlim(sac.b, sac.e)
        ax1.set_xlim(sac.b, sac.e)
        ax1.set_ylim(y.min(), y.max())
        ax2.set_ylim(y.min(), y.max())
        ax3.set_ylim(y.min(), y.max())
        plt.axis('tight')

        def onselect1(xmin, xmax):
            """
            Redraw subplots 2 and 3 after selecting
            window in subplot 1.
            """
            indmin, indmax = npy.searchsorted(x, (xmin, xmax))
            indmax = min(len(x)-1, indmax)
            thisx = x[indmin:indmax]
            thisy = y[indmin:indmax]
            line2.set_data(thisx, thisy)
            line3.set_data(thisx, thisy)
            ax2.set_xlim(thisx[0], thisx[-1])
            ax3.set_xlim(thisx[0], thisx[-1])
            ymin = thisy.min()
            ymax = thisy.max()
            ax2.set_ylim(ymin, ymax)
            ax3.set_ylim(ymin, ymax)
            fig.canvas.draw()

        def onselect2(lxmin, lxmax):
            """
            Redraw subplot 3 after selecting
            windown in subplot 2.
            """
            global xmin, xmax   # I do not know how to get return values other than like this
            xmin = lxmin
            xmax = lxmax
            indmin, indmax = npy.searchsorted(x, (xmin, xmax))
            indmax = min(len(x)-1, indmax)
            thisx = x[indmin:indmax]
            thisy = y[indmin:indmax]
            line3.set_data(thisx, thisy)
            ax3.set_xlim(thisx[0], thisx[-1])
            ax3.set_ylim(thisy.min(), thisy.max())
            fig.canvas.draw()

        span1 = SpanSelector(ax1, onselect1,  'horizontal', useblit=True,
                            rectprops=dict(alpha=0.5, facecolor='red'))
        span2 = SpanSelector(ax2, onselect2, 'horizontal', useblit=True,
                            rectprops=dict(alpha=0.5, facecolor='yellow'))
        fig.canvas.draw()
        plt.show()

    class Index():
        ind = 0
        def save(self, event):
            sac.__setattr__(t1, xmin)
            sac.__setattr__(t2, xmax)
            print('{0}: saved {1}={2}, {3}={4}'.format(sac.fh.name, t1, xmin, t2, xmax))

        def next(self, event):
            self.ind += 1
            try:
                infile = files[self.ind]
            except IndexError:
                self.ind = 0
                infile = files[self.ind]
            sac, x, y = readsac(infile)
            plotData(sac, x, y,  args)

        def prev(self, event):
            self.ind -= 1
            try:
                infile = files[self.ind]
            except IndexError:
                self.ind = -1
                infile = files[self.ind]
            sac, x, y = readsac(infile)
            plotData(sac, x, y,  args)

        def quit(self, event):
            self.close()

    # read params and create subplots
    t1, t2 = args.header
    files = args.sacfile
    fig = plt.figure(figsize=(10,8))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    plt.subplots_adjust(bottom=0.2)

    callback = Index()

    # Create Buttons
    axprev = plt.axes([0.59, 0.05, 0.1, 0.075])
    axnext = plt.axes([0.70, 0.05, 0.1, 0.075])
    axquit = plt.axes([0.81, 0.05, 0.1, 0.075])
    axsave = plt.axes([0.48, 0.05, 0.1, 0.075])
    bnprev = Button(axprev, 'Prev')
    bnprev.on_clicked(callback.prev)
    bnnext = Button(axnext, 'Next')
    bnnext.on_clicked(callback.next)
    bnquit = Button(axquit, 'Quit')
    bnquit.on_clicked(callback.quit)
    bnsave = Button(axsave, 'Save')
    bnsave.on_clicked(callback.save)

    sac, x, y = readsac(files[0])
    plotData(sac, x, y, args)

def main():
    args = getargs()
    StartPlot(args)

if __name__ == "__main__":
    main()
