import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import csv
import numpy as np

def genPlot(arg):
    if  arg=='zad1':
        gs = gridspec.GridSpec(2, 1)
        x = []
        y = []
        with open('f-1.csv', 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter=',')
            for row in plots:
                x.append(float(row[0]))
                y.append(float(row[1]))
        ax = pl.subplot(gs[0, 0])  # row 0, col 1
        pl.title('Sygnal prosty')
        pl.ylabel('s(t)')
        pl.xlabel('Czas[s]')
        pl.plot(x, y, 'o', color='black', markersize=0.5)

        x = []
        y = []
        with open('r-1.csv', 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter=',')
            for row in plots:
                x.append(float(row[0]))
                y.append(float(row[1]))
        ax = pl.subplot(gs[1, 0])  # row 0, col 1
        pl.title('Sygnal prosty po DFT / IDFT')
        pl.ylabel('s(t)')
        pl.xlabel('Czas[s]')
        pl.plot(x, y, 'o', color='black', markersize=0.5)
        pl.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.5, hspace=0.5)
        pl.savefig('Zad1-4.png', dpi=500, bbox_inches='tight', pad_inches=0.1)
    else:
        x = []
        y = []
        with open(arg + '-1.csv', 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter=',')
            for row in plots:
                x.append(float(row[0]))
                y.append(float(row[1]))
        gs = gridspec.GridSpec(2, 2)

        pl.figure()
        ax = pl.subplot(gs[0, 0]) # row 0, col 0
        pl.title('Sygnal wejsciowy')
        pl.ylabel(arg + '(t)')
        pl.xlabel('Czas[s]')
        #pl.plot(x, y, '--ok', color='black', markersize=1, linewidth=0.5);
        pl.plot(x, y, 'o', color='black', markersize=0.5)

        x = []
        y = []
        with open(arg + '-2.csv', 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter=',')
            for row in plots:
                x.append(float(row[0]))
                y.append(float(row[1]))
        ax = pl.subplot(gs[0, 1]) # row 0, col 1
        pl.title('Widmo amplitudowe')
        pl.ylabel('A[V]')
        pl.xlabel('f[Hz]')
        #pl.plot(x, y)
        pl.stem(x, y, use_line_collection="True", markerfmt=" ")

        x = []
        y = []
        with open(arg + '-3.csv', 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter=',')
            for row in plots:
                x.append(float(row[0]))
                y.append(float(row[1]))
        ax = pl.subplot(gs[1, :]) # row 1, span all columns
        pl.title('Widmo amplitudowe w skali decybelowej')
        pl.ylabel('A[dB]')
        pl.xlabel('f[Hz]')
        pl.stem(x, y, use_line_collection="True", markerfmt=" ")
        pl.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.5, hspace=0.5)
        pl.savefig(arg + '.png', dpi=500,  bbox_inches='tight', pad_inches=0.1)

    pl.cla()
    print("Zakonczono " + arg + ".")

genPlot('zad1')
genPlot('s')
genPlot('x')
genPlot('y')
genPlot('z')
genPlot('u')
genPlot('v')
genPlot('p1')
genPlot('p2')
genPlot('p3')