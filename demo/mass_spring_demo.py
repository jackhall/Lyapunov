#!/usr/bin/python

import time
import lyapunov
import solvers
import numpy
import matplotlib.pyplot as plt


class MassSpringDemo(object):
	"""
	Mass spring damper system.
	k = b = m = 1.0
	No disturbances or control.
	"""
	def __init__(self):
		self.state = (1.0, 1.0) 
		self.u = lambda : 0.0

	def __call__(self):
		x, v = self.state
		return (v, -v - x + self.u())

	def plot(self):
		plt.figure()
		try:
			plt.plot(self.x_out[:,0], self.x_out[:,1])
			plt.show()
		except AttributeError:
			plt.close()


sys = MassSpringDemo()
t_in = [0.1, 0.2]
sys.time = 0.0
stepper = solvers.runge_kutta4(sys, t_in)
print "No Events - mass spring damper system"
print "Step 0:"
print "state ", sys.state, "time ", sys.time
print "slope ", sys()

stepper.step(0.1)
print "Step 1:"
print "state ", sys.state, "time ", sys.time
print "slope ", sys()

stepper.step(0.2)
print "Step 2:"
print "state ", sys.state, "time ", sys.time
print "slope ", sys()


class SubsystemDemo(object):
	"""A mass-spring-damper controlled by a PID."""
	def __init__(self):
		self.plant = MassSpringDemo()
		self.control = lyapunov.PID(Ki=1)
		self.control.y = lambda : self.plant.state
		self.control.r = self.reference
		self.plant.u = self.control.u
		self.time = 0.0

	@property
	def state(self):
		return self.control.state + self.plant.state

	@state.setter
	def state(self, x):
		self.control.state = (x[1],)
		self.plant.state = x[1:]

	def __call__(self):
		return self.control() + self.plant()

	def reference(self):
		return 0.0 if self.time < 2.0 else 1.0

	def plot(self):
		self.plant.t_out = self.t_out
		self.plant.x_out = self.x_out[:,1:]
		self.plant.plot()


sys4 = SubsystemDemo()
t_in = numpy.linspace(0.0, 3.0, 200)
print "\nMass-Spring-Damper w/PID control"
print "initial state", sys4.state
start = time.clock()
record = lyapunov.simulate(sys4, t_in)
sys4.x_out, sys4.t_out = numpy.array(record.state), numpy.array(record.time)
print "time elapsed", time.clock() - start
sys4.plot()
