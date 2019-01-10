#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np              # Numerical library
from scipy import *             # Load the scipy functions
from control.matlab import *    # Load the controls systems library
from matplotlib import pyplot as plt

def drange(begin, end, step):
    n = begin
    end = end + step
    while n < end:
     yield n
     n += step


foot = [[0, 0, 0], [0.6, 0.1, 0.06], [0.9, 0.2, -0.06], [1.2, 0.3, 0.06], [1.5, 0.4, -0.06],[1.8, 0.5, 0.06], [2.4, 0.6, -0.06], [3.0, 0.7, 0.0],[100,0,0]]
forward_period = 1.0
calculate_period = 4.0

dt = 0.01
zh = 0.27
g  = 9.8

A = np.matrix([[0, 1, 0],
               [0, 0, 1],
               [0, 0, 0]])

B = np.matrix([[0],
               [0],
               [1]])

C = np.matrix([[1, 0, (-zh/g)]])

D = 0

sys = ss(A,B,C,D)
sys_d = c2d(sys, dt)
A_d, B_d, C_d, D_d = ssdata(sys_d)

E_d  = np.matrix([[dt],[1], [0]])
Zero = np.matrix([[0], [0], [0]])
Phai = np.r_[np.c_[1, -C_d * A_d], np.c_[Zero, A_d]]
G    = np.r_[-C_d*B_d, B_d]
GR   = np.matrix([[1], [0], [0], [0]])
Gd   = np.r_[-C_d*E_d, E_d]
Q = np.zeros((4,4))
Q[0, 0] = pow(10,8)
H = 1

P = dare(Phai, G, Q, H)[0]
F = pow(-(H+(G.transpose())*P*G),-1)*G.transpose()*P*Phai

x = np.matrix([[0],[0],[0]])
y = np.matrix([[0],[0],[0]])
xp = x;
yp = y;

t = np.arange(0,calculate_period + dt,dt)
i = 1;
n = 0;
prefx = []
prefy = []

x0 = []
y0 = []
x1 = []
y1 = []
x2 = []
y2 = []
times = []

for tt in drange(0, calculate_period+forward_period+1-dt, dt):
    time = float(format(round(tt, 3)))
    if(time == foot[n][0]):
        prefx.append(foot[n][1])
        prefy.append(foot[n][2])
        #print(i, time, foot[n][1])
        n = n + 1
    else:
        prefx.append(prefx[i-2])
        prefy.append(prefy[i-2])
    i = i + 1

i = 0
ux,uy = 0.0, 0.0
xi = (np.eye(4)-G*pow(H+G.transpose()*P*G,-1)*G.transpose()*P)*Phai

for ttn in drange(0,calculate_period-dt,dt):
    tt = float(format(round(ttn, 3)))
    px = C_d*x
    py = C_d*y
    ex = prefx[i] - px
    ey = prefy[i] - py
    X = np.r_[ex,x - xp]
    Y = np.r_[ey,y - yp]
    xp = x
    yp = y
    dux = F * X
    j = 0

    for tttn in drange(tt,tt+ forward_period-dt, dt):
        j = j + 1
        if (prefx[i+j] - prefx[i+j-1]) != 0:
            f  = -pow((H+G.transpose()*P*G),-1)*G.transpose()*pow(xi.transpose(),j-1)*P*GR
            dux = dux + f * (prefx[i+j] - prefx[i+j-1])
    ux += dux
    duy = F * Y
    j = 0

    for tttn in drange(tt,tt+ forward_period-dt, dt):
        j = j + 1
        if (prefy[i+j] - prefy[i+j-1]) != 0:
            f  = -pow((H+G.transpose()*P*G),-1)*G.transpose()*pow(xi.transpose(),j-1)*P*GR
            duy = duy + f * (prefy[i+j] - prefy[i+j-1])
    uy += duy

    dx = 0
    dy = 0
    x = A_d  * x + B_d * ux + E_d * dx * dt
    y = A_d * y + B_d * uy + E_d * dy * dt;
    x0.append(float(x[0, 0]))
    y0.append(float(y[0, 0]))
    x1.append(float(prefx[i]))
    y1.append(float(prefy[i]))
    x2.append(float(px))
    y2.append(float(py))
    times.append(tt)
    i = i + 1

g = plt.subplot(2,1,1)
g.plot(x0, y0, color="red",label="$COM$")
g.plot(x1, y1, color="green",label="$refZMP$")
g.plot(x2, y2, "*",label="$ZMP$")
g.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=8)
r = plt.subplot(2,1,2)

r.plot(times, x0, color="red",label="$COM X$")
r.plot(times, x1, color="blue",label="$refZMPX$")
r.plot(times, x2, color="lime",label="$ZMP X$")
r.plot(times, y0, color="green",label="$COM Y$")
r.plot(times, y1, color="cyan",label="$refZMPY$")
r.plot(times, y2, color="violet",label="$ZMP Y$")
r.legend(bbox_to_anchor=(1, 0.9), loc='upper right', borderaxespad=0, fontsize=8)
plt.savefig('preview_control.png')
plt.show()
