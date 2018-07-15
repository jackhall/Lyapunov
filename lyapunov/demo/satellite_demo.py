#!/usr/bin/python

#Lyapunov: a library for integrating nonlinear dynamical systems
#Copyright (C) 2013-2018 John Wendell Hall
#
#The author may be reached at jackwhall7@gmail.com.

import time
import lyapunov
import numpy
import matplotlib.pyplot as plt


class SatelliteDemo(object):
	"""Double integrator with linear switching mode."""
	def __init__(self):
		self.state = 0.0, (1.0, 1.0)
		self.u = lambda: (1 if self.s() < 0 else -1)

	state = lyapunov.state_property(xname="_state")

	def s(self):
		"""event function"""
		x, v = self.state.x
		return v + 0.5*x

	def u_margin(self):
		"""event function: negative when control limits exceeded"""
		return 1.0 - abs(self.u_effective())

	def u_effective(self):
		x, v = self.state.x
		return v**2 / x

	def __call__(self):
		_, v = self.state.x
		return (v, self.u())


def plot_satellite(record):
    plt.figure()
    x = numpy.array(record.x)
    xmin, xmax = min(x[:,0]), max(x[:,0])
    ymin, ymax = min(x[:,1]), max(x[:,1])
    x_mode = numpy.array([xmin, xmax])
    lam = numpy.array([-0.5*xmin, -0.5*xmax])
    plt.plot(x_mode, lam)
    plt.plot(x[:,0], x[:,1])
    plt.show()


def run_satellite_demo_noevents():
    #Brute force could work well, but with a small step size...
    system = SatelliteDemo()
    #sol = lyapunov.Solver(system)
    t_in = numpy.linspace(0.0, 3.0, 31)
    record = lyapunov.Recorder(system)
    stepper = lyapunov.runge_kutta4(system, t_in)
    print("\nNo Events, No Sliding - satellite control")
    print("initial state", system.state)
    start = time.clock()
    #record = sol.simulate(t_in)
    for t, events in stepper:
        record.log()
    print("time elapsed", time.clock() - start)
    plot_satellite(record)


def run_satellite_demo_events():
    #With events and sliding is more accurate, and faster than brute force.
    system = SatelliteDemo()
    system.events = [system.s]
    system.u = lambda: -1 #replaces u_naive
    t_in = numpy.linspace(0.0, 3.0, 31)
    record = lyapunov.Recorder(system)
    stepper = lyapunov.runge_kutta4(system, t_in)
    print("\nWith Events and Sliding - satellite control")
    start = time.clock()
    for t, events in stepper:
        record.log(events)
        if len(events) > 0:
            if system.s in events and system.u_margin() > 0:
                system.u = system.u_effective
                system.events = [system.u_margin]
            else:
                stepper.step_through()
                system.u = (lambda: 1) if system.u_effective() > 1 else (lambda: -1)
                system.events = [system.s]
    print("time elapsed ", time.clock() - start)
    plot_satellite(record)


if __name__ == "__main__":
    run_satellite_demo_events()
