#!/usr/bin/python

import lyapunov
import numpy
from scipy.interpolate import interp1d

pendubot_data_begin = {"pendubot_run1.dat": 210, 
              "pendubot_run2.dat": 230,
              "pendubot_run3.dat": 170}

class PendubotReader(object):
    def __init__(self, filename):
        data = numpy.loadtxt(filename, 
                             skiprows=pendubot_data_begin[filename], 
                             usecols=range(5))
        rows = data.shape[0]
        self.time_data = numpy.linspace(0, rows*0.01, rows)
        self.theta1 = interp1d(time_data, data[:,1], kind=cubic)
        self.theta2 = interp1d(time_data, data[:,2], kind=cubic)
        #self.omega1 = interp1d(data[:,0], data[:,3], kind=cubic)
        #self.omega2 = interp1d(data[:,0], data[:,4], kind=cubic)
        self.state = 0.0, ()

    state = lyapunov.state_property(tname='time')

    def theta(self):
        return self.theta1(self.time), self.theta2(self.time)


class PendubotObserver(object):
    def __init__(self):
        self.state = 0.0, (0.0,)*6
        self.events = 
        self.y
       
    state = lyapunov.state_property()

    def __call__(self):
        pass

def run_pendubot_demo():
    filename = "pendubot_run1.dat"
    reader = PendubotReader(filename)
    observer = PendubotObserver(filename)
    pendubot = lyapunov.ParallelSystems([reader, observer])
    observer.y = reader.theta
    observer.state = 0.0, reader.theta + (0.0, 0.0) + (0.1, 0.1) 

    stepper = lyapunov.adams_bashforth4(pendubot, reader.time_data)
    for t, events in stepper:
        if events:
            
    

