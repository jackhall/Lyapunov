#!/usr/bin/python

#Lyapunov: a library for integrating nonlinear dynamical systems
#Copyright (C) 2013  John Wendell Hall
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#The author may be reached at jackhall@utexas.edu.

import time
import lyapunov
import numpy
import matplotlib.pyplot as plt


def run_filter_demo():
    fil = lyapunov.Filter((1.0, 3.0, 3.0))
    ref = lyapunov.StepSignal()
    fil.signal = lambda : ref.value
    fil.state = 0.0, (0.0,)*3
    system = lyapunov.ParallelSystems([ref, fil])
    record = lyapunov.Recorder(system)
    stepper = lyapunov.dormand_prince(system, 10.0)
    print "\nFilter"
    print "initial state", system.state
    start = time.clock()
    for t, events in stepper:
        record.log(events)
        if events:
            stepper.step_across()
            ref.update()
    print "time elapsed", time.clock() - start
    plt.figure()
    try:
        plt.plot(record.t, numpy.array(record.x)[:,0])
        plt.show()
    except:
        plt.close()


if __name__ == "__main__":
    run_filter_demo()

