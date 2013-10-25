#!/usr/bin/python

import math
import numpy
import sympy as sym
from scipy.interpolate import interp1d
from functools import partial
import lyapunov

theta1, theta2, h1, h2 = sym.symbols("theta1 theta2 h1 h2")
m1, m2, J1, J2 = sym.symbols("m1 m2 J1 J2")
l1, l2, L1 = sym.symbols("l1 l2 L1")
b1, b2 = sym.symbols("b1 b2")
g = sym.symbols("g")
theta1dot = h1 / (m1*l1**2 + J1)
theta2dot = (h2 - (m2*L1*(L1 + l2*sym.cos(theta2)) + J2)*theta1dot) / (J2 
            + m2*l2*sym.cos(theta2)*(L1 + l2*sym.cos(theta2)))
h1dot = m1*g*l1*sym.sin(theta1) + m2*g*(L1*sym.sin(theta1) 
        + l2*sym.sin(theta1+theta2)) - b1*theta1dot
h2dot = m2*g*l2*sym.sin(theta1+theta2) - b2*theta2dot
x = sym.Matrix([theta1, theta2, h1, h2, b1, b2])
f = sym.Matrix([theta1dot, theta2dot, h1dot, h2dot, sym.S(0), sym.S(0)])
h = sym.Matrix([theta1, theta2])

df = f.jacobian(x)
dh = h.jacobian(x)
L11, L12, L13, L14, L15, L16 = sym.symbols("L11 L12 L13 L14 L15 L16")
L21, L22, L23, L24, L25, L26 = sym.symbols("L21 L22 L23 L24 L25 L26")
L = sym.Matrix([[L11, L21],[L12, L22],[L13, L23],[L14, L24],[L15, L25],[L16, L26]])

char_eqn = sym.Matrix.charpoly(df - L*dh, lambda expr: sym.collect(sym.expand(expr)))

class Reader(object):
    def __init__(self, filename):
        data = numpy.loadtxt(filename, 
                             usecols=range(5))
        rows = data.shape[0]
        self.time_data = numpy.linspace(0, rows*0.01, rows)
        theta1_data = numpy.unwrap(data[:,1])
        theta2_data = numpy.unwrap(data[:,2])
        self.theta1 = interp1d(time_data, theta1_data, kind=cubic)
        self.theta2 = interp1d(time_data, theta2_data, kind=cubic)
        #self.omega1 = interp1d(data[:,0], data[:,3], kind=cubic)
        #self.omega2 = interp1d(data[:,0], data[:,4], kind=cubic)
        self.state = 0.0, ()

    state = lyapunov.state_property(tname='time')

    def theta(self):
        return self.theta1(self.time), self.theta2(self.time)


class Observer(object):
    def __init__(self):
        self.state = 0.0, (0.0,)*6
        self.y = None
       
    state = lyapunov.state_property()

    def __call__(self):
        #return (ydot - output_error*L)
        pass


def run_pendubot_demo():
    filename = "pendubot_run1.dat"
    reader = Reader(filename)
    observer = Observer(filename)
    pendubot = lyapunov.ParallelSystems([reader, observer])
    observer.y = reader.theta
    observer.state = 0.0, reader.theta + (0.0, 0.0) + (0.1, 0.1) 

    record = lyapunov.Recorder(pendubot)
    stepper = lyapunov.adams_bashforth4(pendubot, reader.time_data)
    for t, _ in stepper:
        record.log()                        
                   
    print "b1 ~=", observer.state.x[-2]
    print "b2 ~=", observer.state.x[-1]

