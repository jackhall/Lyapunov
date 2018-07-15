#!/usr/bin/python

import math
import numpy, numpy.linalg
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
        statelist = [theta1, theta2, h1, h2, b1, b2]
        theta1dot = h1 / (m1*l1**2 + J1)
        theta2dot = (h2 - (m2*L1*(L1 + l2*sym.cos(theta2)) + J2)*theta1dot) / (J2
                    + m2*l2*sym.cos(theta2)*(L1 + l2*sym.cos(theta2)))
        h1dot = m1*g*l1*sym.sin(theta1) + m2*g*(L1*sym.sin(theta1)
                + l2*sym.sin(theta1+theta2)) - b1*theta1dot
        h2dot = m2*g*l2*sym.sin(theta1+theta2) - b2*theta2dot
        self.xsym = sym.Matrix(statelist)
        o = sym.S(0)
        f = sym.Matrix([theta1dot, theta2dot, h1dot, h2dot, o, o]).subs(parameters)
        df = f.jacobian(self.xsym)
        self.f = sym.lambdify(statelist, f)
        self.df = sym.lambdify(statelist, df)

        #set poles
        A = sym.Matrix([[a,sym.S(-1),j,o,o,o],
                        [b,sym.S(-1),k,l,o,o],
                        [c,sym.S(-1),m,o,n,o],
                        [d,sym.S(-1),p,q,o,r],
                        [e,sym.S(-1),o,o,o,o],
                        [i,sym.S(-1),o,o,o,o]])
        self.desired_coefficients = numpy.array([[729], [1458], [1215], [540], [135], [18]])
        char_poly = (sym.symbols("_lambda")*sym.eye(6) - A).det_bareis()
        char_poly = [char_poly.coeff("_lambda", index) for index in range(6)]
        char_poly = sym.Matrix([[sym.diff(row, var) for var in A[:,0]]
                                                    for row in char_poly])
        char_poly = sym.lambdify((j,k,l,m,n,p,q,r), char_poly)
        self.char_poly = lambda x: numpy.matrix(char_poly(*x))

    def __call__(self):
        dfeval = self.df(*self.state.x)
        #bind current values of df to charateristic polynomial
        gradientsubs = (dfeval[0,2], dfeval[1,2], dfeval[1,3], dfeval[2,2],
                        dfeval[2,4], dfeval[3,2], dfeval[3,3], dfeval[3,5])
        char_poly = self.char_poly(gradientsubs)
        if numpy.linalg.cond(char_poly) < 1000:
            col1 = numpy.linalg.solve(char_poly, self.desired_coefficients)
            L = [[-col1[0], 1],
                 [-col1[1], 1 - dfeval[1,1]],
                 [-col1[2] - dfeval[2,0], 1 - dfeval[2,1]],
                 [-col1[3] - dfeval[3,0], 1 - dfeval[3,1]],
                 [-col1[4], 1],
                 [-col1[5], 1]]
        else:
            L = [[0,0] for index in range(6)]

        #compute error of current output estimate given actual outputs
        error = [y - yhat for y, yhat in zip(self.y(), self.state.x[:2])]
        return tuple(float(xdoti - Lrow[0]*error[0] - Lrow[1]*error[1])
                     for xdoti, Lrow in zip(self.f(*self.state.x), L))


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

    print("b1 ~=", observer.state.x[-2])
    print("b2 ~=", observer.state.x[-1])
    return record
