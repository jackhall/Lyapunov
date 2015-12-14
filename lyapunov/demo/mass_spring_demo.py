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


class MassSpringDemo(object):
    """
    Mass spring damper system.
    k = b = m = 1.0
    No disturbances or control.
    """
    def __init__(self):
        self.state = 0.0, (1.0, 1.0)
        self.u = lambda : 0.0

    state = lyapunov.state_property(xname="_state")

    def __call__(self):
        x, v = self._state
        return (v, -v - x + self.u())


def print_two_steps(system = MassSpringDemo()):
    """ System must make use of lyapunov.State. """
    t_in = [0.1, 0.2]
    stepper = lyapunov.euler(system, t_in)
    print "No Events - mass spring damper systemtem"
    print "Step 0:"
    print "time", system.state.t, "| state ", system.state.x
    print "slope", system()

    stepper.next()
    print "Step 1:"
    print "time", system.state.t, "| state ", system.state.x
    print "slope", system()

    stepper.next()
    print "Step 2:"
    print "time", system.state.t, "| state ", system.state.x
    print "slope", system()


class SubsystemDemo(lyapunov.ParallelSystems):
    """A mass-spring-damper controlled by a PID."""
    def __init__(self):
        self.plant = MassSpringDemo()
        self.control = lyapunov.PID(Ki=1)
        self.reference = lyapunov.StepSignal(step_time=4.0)
        self.control.y = lambda: self.plant.state[1]
        self.control.r = lambda: self.reference.value
        self.plant.u = self.control.u
        lyapunov.ParallelSystems.__init__(self, [self.reference, 
                                                 self.control, 
                                                 self.plant])


def run_subsystem_demo():
    system = SubsystemDemo()
    record = lyapunov.Recorder(system)
    stepper = lyapunov.adams_bashforth3(system, numpy.linspace(0.0, 8.0, 100))
    #stepper = lyapunov.cash_karp(system, 8.0)

    print "\nMass-Spring-Damper w/PID control"
    print "initial state", system.state
    start = time.clock()
    count = 0
    for t, events in stepper:
        if events:
            stepper.step_across()
            system.reference.update()
        record.log(events)
    print "time elapsed", time.clock() - start

    x_out = numpy.array(record.x)
    plt.figure()
    plt.plot(x_out[:,0], x_out[:,1])
    plt.show()


if __name__ == "__main__":
    run_subsystem_demo();

