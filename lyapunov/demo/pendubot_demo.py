#!/usr/bin/python

import math
import numpy
import sympy as sym
from sympy.abc import a,b,c,d,e,g,i,j,k,l,m,n,p,q,r
from scipy.interpolate import interp1d
from functools import partial
import lyapunov

class Reader(object):
    def __init__(self, filename):
        data = numpy.loadtxt(filename, 
                             usecols=range(5))
        rows = data.shape[0]
        self.time_data = numpy.linspace(0, rows*0.01, rows)
        theta1_data = numpy.unwrap(data[:,1])
        theta2_data = numpy.unwrap(data[:,2])
        self.theta1 = interp1d(self.time_data, theta1_data, kind='cubic')
        self.theta2 = interp1d(self.time_data, theta2_data, kind='cubic')
        #self.omega1 = interp1d(data[:,0], data[:,3], kind=cubic)
        #self.omega2 = interp1d(data[:,0], data[:,4], kind=cubic)
        self.state = 0.0, ()
    state = lyapunov.state_property(tname='time')
    def theta(self):
        return self.theta1(self.time), self.theta2(self.time)


class Observer(object):
    state = lyapunov.state_property()
    def __init__(self, y):
        self.state = 0.0, (0.0,)*6
        self.y = y #function returning actual angles
        #define pendubot model symbolically
        m1, m2, J1, J2 = sym.symbols("m1 m2 J1 J2")
        l1, l2, L1 = sym.symbols("l1 l2 L1")
        parameters = {g: 9.81,
                      m1: 1,
                      m2: 1,
                      J1: .1,
                      J2: .1,
                      l1: .2,
                      l2: .2,
                      L1: .3}
        theta1, theta2, h1, h2, b1, b2 = sym.symbols("theta1 theta2 h1 h2 b1 b2")
        theta1dot = h1 / (m1*l1**2 + J1)
        theta2dot = (h2 - (m2*L1*(L1 + l2*sym.cos(theta2)) + J2)*theta1dot) / (J2 
                    + m2*l2*sym.cos(theta2)*(L1 + l2*sym.cos(theta2)))
        h1dot = m1*g*l1*sym.sin(theta1) + m2*g*(L1*sym.sin(theta1) 
                + l2*sym.sin(theta1+theta2)) - b1*theta1dot
        h2dot = m2*g*l2*sym.sin(theta1+theta2) - b2*theta2dot
        self.xsym = sym.Matrix([theta1, theta2, h1, h2, b1, b2])
        o = sym.S(0)
        self.f = sym.Matrix([theta1dot, theta2dot, h1dot, h2dot, o, o]).subs(parameters)
        self.df = self.f.jacobian(self.xsym)
        #set poles
        A = sym.Matrix([[a,sym.S(-1),j,o,o,o],
                        [b,sym.S(-1),k,l,o,o],
                        [c,sym.S(-1),m,o,n,o],
                        [d,sym.S(-1),p,q,o,r],
                        [e,sym.S(-1),o,o,o,o],
                        [i,sym.S(-1),o,o,o,o]])
        desired_coeff = [729, 1458, 1215, 540, 135, 18]
        self.char_poly = (sym.symbols("_lambda")*sym.eye(6) - A).det_bareis()
        self.char_poly = [self.char_poly.coeff("_lambda", index) - desired 
                          for index, desired in enumerate(desired_coeff)]
    def __call__(self):
        #bind current values of x to df
        xsubs = {xsym: xval for xsym, xval in zip(self.xsym, self.state.x)}
        dfeval = sym.Matrix(self.df).subs(xsubs)
        #bind current values of df to charateristic polynomial
        gradsubs = {j:dfeval[0,2], k:dfeval[1,2], l:dfeval[1,3], m:dfeval[2,2], 
                    n:dfeval[2,4], p:dfeval[3,2], q:dfeval[3,3], r:dfeval[3,5]}
        char_poly_desired = map(lambda term: term.subs(gradsubs), list(self.char_poly))
        #solve for observer gains
        col1 = sym.solve(char_poly_desired, [a,b,c,d,e,i])
        if not col1:
            #observability lost
            pass
        L = [[-col1[a], 1],
             [-col1[b], 1 - dfeval[1,1]],
             [-col1[c] - dfeval[2,0], 1 - dfeval[2,1]],
             [-col1[d] - dfeval[3,0], 1 - dfeval[3,1]],
             [-col1[e], 1],
             [-col1[i], 1]]
        #bind xsubs to self.f
        xhatdot = sym.Matrix(self.f).subs(xsubs)
        #compute error of current output estimate given actual outputs
        error = [y - yhat for y, yhat in zip(self.y(), self.state.x[:2])]
        return tuple(float(xdoti.evalf() - Lrow[0]*error[0] - Lrow[1]*error[1])
                     for xdoti, Lrow in zip(xhatdot, L))


def run_pendubot_demo():
    filename = "pendubot_run1.dat"
    reader = Reader(filename)
    observer = Observer(reader.theta)
    pendubot = lyapunov.ParallelSystems([reader, observer])
    observer.y = reader.theta
    observer.state = 0.0, reader.theta() + (0.1, 0.1) + (0.1, 0.1) 

    record = lyapunov.Recorder(pendubot)
    stepper = lyapunov.adams_bashforth4(pendubot, reader.time_data)
    for t, _ in stepper:
        record.log()                        
                  
    print "b1 ~=", observer.state.x[-2]
    print "b2 ~=", observer.state.x[-1]
    
    record.plot()

