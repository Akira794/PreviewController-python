#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np              # Numerical library
from scipy import *             # Load the scipy functions
from control.matlab import *    # Load the controls systems library
from matplotlib import pyplot as plt

class CurveFitting():

    def __init__(self, _x, _t, _dim):
        self.x = _x
        self.t = _t
        self.dim = _dim
        self.W = 0

    def calc_polynomial_fitting(self):
        A = np.zeros((self.dim+1,self.dim+1))
        i = 0
        for i in range(self.dim+1):
            for j in range(self.dim+1):
                temp = np.power(self.x,i+j).sum()
                if( i == j):
                    temp += self.dim
                A[i, j] = temp

        T = []
        for n in range(self.dim+1):
            T.append(np.dot(np.power(self.x,n),self.t).sum())

        #print(A)
        #print(T)
        self.W = np.linalg.solve(A,T)
        #print(self.W)

    def out_y(self, x_range):
        result = self.W[0]
        for k in range(self.dim+1):
            result += self.W[k]*pow(x_range,k)
        return result

x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
y = [0.0, 1.0, 2.0, 4.5, 8.0, 12.5, 20.0, 30.5, 35, 40.5, 50]
CF = CurveFitting(x,y,5)
CF.calc_polynomial_fitting()
print(CF.W)
p = []
for i in x:
    p.append(CF.out_y(i))

plt.plot(x,y,color="red",label="$Base$")
plt.plot(x,p,color="green",label="$CurveFitting$")
plt.show()
