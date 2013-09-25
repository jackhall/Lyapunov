#!/usr/bin/python

import math
import numpy
from scipy.interpolate import interp1d
import lyapunov

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

