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

dist = [0.6, 0.05, ( 3.14/2.0) ]
max_step_x = 0.1
max_step_y = 0.05
max_step_w = 0.2;
period = 0.3
foot_y = 0.06

step_num_x = int(abs(dist[0])/max_step_x)
step_num_y = int(abs(dist[1])/max_step_y)
step_num_w = int(abs(dist[2])/max_step_w)
step_num = max([step_num_x, step_num_y, step_num_w])

if (step_num_x == step_num):
    step_x = sign(dist[0])*min(abs(dist[0]), max_step_x)
    step_y = dist[1]/abs(dist[0])*max_step_x
    step_w = dist[2]/abs(dist[0])*max_step_x
elif (step_num_y == step_num):
    step_x = dist[0]/abs(dist[1])*max_step_y
    step_y = sign(dist[1])*min(abs(dist[1]), max_step_y)
    step_w = dist[2]/abs(dist[1])*max_step_y
elif (step_num_w == step_num):
    step_x = dist[0]/abs(dist[2])*max_step_w
    step_y = dist[1]/abs(dist[2])*max_step_w
    step_w = sign(dist[2])*min(abs(dist[2]), max_step_w)

step = np.array([0.0, 0.0, 0.0])
foot = [step]
foot.append(np.array([period, 0.0, foot_y]))
rot = 0

for i in drange(1, step_num, 1):
    shift_y = foot_y * ((divmod(i,2)[1])*2-1)*2*(-1)
    foot_rd = [step_x, step_y+shift_y]
    foot_fd = [foot_rd[0]*cos(rot)-foot_rd[1]*sin(rot), foot_rd[0]*sin(rot)+foot_rd[1]*cos(rot)]
    rot = rot + step_w
    foot.append(np.array([foot[i][0] + period, foot[i][1] + foot_fd[0], foot[i][2] + foot_fd[1]]))

foot_rd = [dist[0]-step_x*step_num, dist[1]-step_y*step_num-shift_y]
foot_fd = [foot_rd[0]*cos(rot)-foot_rd[1]*sin(rot), foot_rd[0]*sin(rot)+foot_rd[1]*cos(rot)]
foot.append(np.array([foot[i+1][0] + period, foot[i+1][1] + foot_fd[0], foot[i+1][2] + foot_fd[1]]))

rot = dist[2]
foot_rd = [0,shift_y/2]
foot_fd = [foot_rd[0]*cos(rot)-foot_rd[1]*sin(rot), foot_rd[0]*sin(rot)+foot_rd[1]*cos(rot)]
foot.append(np.array([foot[i+2][0] + period, foot[i+2][1] + foot_fd[0], foot[i+2][2] + foot_fd[1]]))
foot.append(np.array([100, 0, 0]))

#foot = [[0, 0, 0], [0.6, 0.1, 0.06], [0.9, 0.2, -0.06], [1.2, 0.3, 0.06], [1.5, 0.4, -0.06],[1.8, 0.5, 0.06], [2.4, 0.6, -0.06], [3.0, 0.7, 0.0],[100,0,0]]
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
    if(abs(time-foot[n][0])<(dt/2)):
        prefx.append(foot[n][1])
        prefy.append(foot[n][2])
        #print(i, time, foot[n][1])
        n = n + 1
    else:
        prefx.append(prefx[i-2])
        prefy.append(prefy[i-2])
    i = i + 1
"""
for num in range(len(prefx)):
    print(prefx[num], prefy[num])
"""
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
"""
g = plt.subplot(2,1,1)
g.plot(x0, y0, color="red",label="$COM$")
g.plot(x1, y1, color="green",label="$refZMP$")
g.plot(x2, y2, "*",label="$ZMP$")
g.set_aspect('equal')
r = plt.subplot(2,1,2)
r.plot(times, x0, color="red",label="$COM X$")
r.plot(times, x1, color="blue",label="$refZMPX$")
r.plot(times, x2, color="lime",label="$ZMP X$")
r.plot(times, y0, color="green",label="$COM Y$")
r.plot(times, y1, color="cyan",label="$refZMPY$")
r.plot(times, y2, color="violet",label="$ZMP Y$")
plt.savefig('preview_control_rot.png')
plt.show()
"""
step_times = []
step1_x = []
step1_y = []
coefficient_x = []
coefficient_y = []
coefficient_t = []

num = 0

one_step = 30#30
start_step = one_step
for n_step in range(len(foot)):
    #if(n_step > 0 and n_step < len(foot)-2 ):
    if(n_step < len(foot)-1 ):
        #print(n_step*(one_step)+1, (n_step+1)*(one_step)+1)
        for num in range(n_step*(one_step)+1, (n_step+1)*(one_step)+1):
            step1_x.append(x0[num])
            step1_y.append(y0[num])
            step_times.append(times[num])

        res_x = np.polyfit(step_times,step1_x, 5)
        res_y = np.polyfit(step_times,step1_y, 5)
        pl_x = np.poly1d(res_x)
        pl_y = np.poly1d(res_y)
        coefficient_x.append(pl_x)
        coefficient_y.append(pl_y)
        coefficient_t.append(np.array(step_times))

        step_times.clear()
        step1_x.clear()
        step1_y.clear()


for n in range(len(coefficient_t)):
    print(list(coefficient_x[n]), list(coefficient_y[n]))


plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.axes().set_aspect('equal')
plt.plot(x0, y0, color="red",label="$COM$")
plt.plot(x1, y1, color="green",label="$refZMP$")
plt.plot(x2, y2, "*",label="$ZMP$")
plt.legend(bbox_to_anchor=(0.0, 1.0), loc='upper left', borderaxespad=0, fontsize=10)
plt.savefig('previewcontroller_walk_rot_COM.png')
plt.show()

plt.xlabel('t(s)')
plt.ylabel('dist(m)')

p = []

for d_num in range(len(coefficient_t)):
    p = coefficient_t[d_num].tolist()
#    plt.plot(p, coefficient_x[d_num](p), color = list_draw_x[d_num])#color="red")
#    plt.plot(p, coefficient_y[d_num](p), color = list_draw_y[d_num])#color="green")
    plt.plot(p, coefficient_x[d_num](p), color ="red")
    plt.plot(p, coefficient_y[d_num](p), color ="darkblue")#color="red")
plt.plot(p, coefficient_y[d_num](p), color="red",label="$POLY COM X$")
plt.plot(p, coefficient_y[d_num](p), color="darkblue",label="$POLY COM Y$")
#plt.plot(x0, y0, color="red",label="$COM$")
#plt.plot(x1, y1, color="green",label="$refZMP$")
#plt.plot(x2, y2, "*",label="$ZMP$")
#plt.legend(bbox_to_anchor=(0.0, 1.0), loc='upper left', borderaxespad=0, fontsize=10)
#plt.savefig('PC_walk_rot_COM.png')
plt.plot(times, x0, linestyle = "dotted",color="green",label="$COM X$")
#plt.plot(times, x1, color="blue",label="$refZMPX$")
#plt.plot(times, x2, color="lime",label="$ZMP X$")
plt.plot(times, y0, linestyle = "dotted",color="orangered",label="$COM Y$")
#plt.plot(times, y1, color="cyan",label="$refZMPY$")
#plt.plot(times, y2, color="violet",label="$ZMP Y$")
plt.legend(bbox_to_anchor=(0.0, 1.0), loc='upper left', borderaxespad=0, fontsize=10)
#plt.savefig('polynomial_COM_trajectory_txty.png')
plt.show()
