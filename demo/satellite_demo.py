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
import solvers
import numpy
import matplotlib.pyplot as plt


class SatelliteDemo(object):
	"""Double integrator with linear switching mode."""
	def __init__(self):
		self.state = (1.0, 1.0)
		self.u = lambda: (1 if self.s() < 0 else -1)

	def s(self):
		"""event function"""
		x, v = self.state
		return v + 0.5*x
	
	def u_margin(self):
		"""event function: negative when control limits exceeded"""
		return 1.0 - abs(self.u_effective())

	def u_effective(self):
		x, v = self.state
		return v**2 / x

	def __call__(self):
		_, v = self.state
		return (v, self.u())

	def plot(self):
		plt.figure()
		try:
			xmin, xmax = min(self.x_out[:,0]), max(self.x_out[:,0])
			ymin, ymax = min(self.x_out[:,1]), max(self.x_out[:,1])
			x = numpy.array([xmin, xmax])
			lam = numpy.array([-0.5*xmin, -0.5*xmax])
			plt.plot(x, lam)
			plt.plot(self.x_out[:,0], self.x_out[:,1])
			plt.show()
		except AttributeError:
			plt.close()


#Brute force could work well, but with a small step size...
sys3 = SatelliteDemo()
#sol = lyapunov.Solver(sys3)
t_in = numpy.linspace(0.0, 3.0, 31)
record = lyapunov.Recorder(sys3)
stepper = solvers.runge_kutta4(sys3, t_in)
print "\nNo Events, No Sliding - satellite control"
print "initial state", sys3.state
start = time.clock()
#record = sol.simulate(t_in)
for t in stepper:
	record.log()
print "time elapsed", time.clock() - start
sys3.x_out, sys3.t_out = numpy.array(record.state), numpy.array(record.time)
sys3.plot()


#With events and sliding is more accurate, and faster than brute force.
sys2 = SatelliteDemo()
sys2.u = lambda: -1 #replaces u_naive
t_in = numpy.linspace(0.0, 3.0, 31)
record = lyapunov.Recorder(sys2)
stepper = solvers.runge_kutta4(sys2, t_in, [sys2.s], 0.0001)
print "\nWith Events and Sliding - satellite control"
start = time.clock()
for t, events in stepper:
	record.log()
	if len(events) > 0:
		if sys2.s in events and sys2.u_margin() > 0:
			sys2.u = sys2.u_effective
			stepper.events = [sys2.u_margin]
		else:
			stepper.step_through()
			sys2.u = (lambda: 1) if sys2.u_effective() > 1 else (lambda: -1)
			stepper.events = [sys2.s]
print "time elapsed ", time.clock() - start
sys2.x_out, sys2.t_out = numpy.array(record.state), numpy.array(record.time)
sys2.plot()

