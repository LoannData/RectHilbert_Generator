#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 19 18:45:03 2019

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.patches as pt
import matplotlib.cm as cm

import hilbertRectangular2D as hi


def hilbertIndexing(U, X, DX, Y, DY) : 
        
#    yimax = int((len(X) - 1)/2.) + int(yidsize/2.)-1
#    yimin = int((len(X) - 1)/2.) - int(yidsize/2.)
    size   = len(U)
    cells  = {}
    loc_index = 0
    for si in range(size) : 
        for xi in range(len(X)) : 
            for yi in range(len(Y)) :
#                if ((X[xi] < U[si][0] and X[xi] + DX > U[si][0]) and 
#                    (Y[yi] < U[si][1] and Y[yi] + DY > U[si][1])) : 
                if (X[xi] == U[si][0] and Y[yi] == U[si][1]) : 
                    temp_cells = {loc_index:[X[xi]-DX/2., Y[yi]-DY/2.]}
                    loc_index += 1
                    cells.update(temp_cells)
    return cells


# Simulation box size 
nx = 4
ny = 6

NX = 2**nx
NY = 2**ny

# X and Y size of the box 
Xmin = 0.
Xmax = 1.
dx = (Xmax - Xmin)/NX


DY   = (Xmax - Xmin)/(float(NX)/float(NY))
Ymax =  DY/2.
Ymin = -DY/2. 



# We define the X and Y vectors
X  = np.linspace(Xmin, Xmax, NX)
Y  = np.linspace(Ymin, Ymax, NY)

dx = X[1] - X[0]
dy = Y[1] - Y[0]

# We get our Hilbert curve 
hilbert = hi.getRHilbert2D(0., 0., NX, NY)


xh = []
yh = []
nhilbert = []

fx = (Xmax - Xmin)/(NX-1)
fy = (Ymax - Ymin)/(NY-1)

for ii in range(len(hilbert)) : 
    xh.append(X[0] + hilbert[ii][0]*fx)
    yh.append(Y[0] + hilbert[ii][1]*fy)
    nhilbert.append([X[0] + hilbert[ii][0]*fx, Y[0] + hilbert[ii][1]*fy])


# Process assignments
nprocess = 16
ncells = NX*NY
sep = int(ncells)/nprocess


cells = hilbertIndexing(nhilbert, X, dx, Y, dy)



# Algorithm to get informations of each cell
cellTable = {}
for ii in range(len(cells)) : 
    Xc = cells.get(ii)[0]
    Yc = cells.get(ii)[1]
    
    for xi in range(len(X)) : 
        for yi in range(len(Y)) : 
            if ((Xc < X[xi] + dx/4. and Xc > X[xi] - dx/4.) and 
                (Yc < Y[yi] + dy/4. and Yc > Y[yi] - dy/4.)) : 
                idX = xi
                idY = yi
    
    tempCellTable = {ii:{"Xc":Xc, "Yc":Yc, "idX":xi, "idY":yi}}
    cellTable.update(tempCellTable)
    
        
        
    
#    temp_cell = {ii:{"Xc":cells.get(ii)[],
#                     "Yc":cells}}








fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111)
#ax.set_xlim(Xmin - 0.2*(Xmax - Xmin), Xmax + 0.2*(Xmax - Xmin))
#ax.set_ylim(0. - 0.2*(Xmax - Xmin), 0. + 0.2*(Xmax - Xmin))

full_box     = pt.Rectangle((X[0]-dx/2.,Y[0]-dy/2.), X[-1] - X[0] + dx, Y[-1] - Y[0] + dy, color='yellow')
ax.add_patch(full_box)

# Plotting the cells borders
for xi in range(len(X)) : 
    for yi in range(len(Y)) :
        ax.plot([X[0]-dx/2., X[-1]+dx/2.],[Y[yi]-dy/2., Y[yi]-dy/2.], color="black", lw=0.5)
    ax.plot([X[0]-dx/2., X[-1]+dx/2.],[Y[-1]+dy/2., Y[-1]+dy/2.], color="black", lw=0.5)
    

for yi in range(len(Y)) : 
    for xi in range(len(X)) : 
        ax.plot([X[xi]-dx/2., X[xi]-dx/2.], [Y[0]-dy/2., Y[-1]+dy/2.], color="black", lw=0.5)
    ax.plot([X[-1]+dx/2., X[-1]+dx/2.], [Y[0]-dy/2., Y[-1]+dy/2.], color="black", lw=0.5)


col_lim = []
for jj in range(nprocess-1) :
    base_limit = jj*sep
    col_limit = (jj+1)*sep
    loc_color = np.random.rand(3,)
    for ii in range(int(base_limit), int(col_limit)) : 
        if (cells.get(ii)) : 
            x0 = cells.get(ii)[0]
            y0 = cells.get(ii)[1]
        loc_cell = pt.Rectangle((x0, y0), dx, dy, color=loc_color)
        ax.add_patch(loc_cell)
#lcol_limit = ncells - col_limit
loc_color = np.random.rand(3,)
for ii in range(int(col_limit), int(ncells)) : 
    if (cells.get(ii)) : 
        x0 = cells.get(ii)[0]
        y0 = cells.get(ii)[1]
    loc_cell = pt.Rectangle((x0, y0), dx, dy, color=loc_color)
    ax.add_patch(loc_cell)


# We plot the Hilbert curve
ax.plot(xh, yh, c="black", lw=2)
ax.set_aspect('equal', 'box')

fig.savefig("./hilbert_example.png")

