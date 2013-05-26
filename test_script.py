#/usr/bin/python

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
import matplotlib.pyplot as plt

sys = lyapunov.SimpleDemo()
sys.time = 0.0
stepper = solvers.Stepper(sys)
print "No Events - mass spring damper system"
print "Step 0:"
print "state ", sys.state, "time ", sys.time
print "error ", stepper.error, "slope ", sys()

stepper.step(0.1)
print "Step 1:"
print "state ", sys.state, "time ", sys.time
print "error ", stepper.error, "slope ", sys()

stepper.step(0.1)
print "Step 2:"
print "state ", sys.state, "time ", sys.time
print "error ", stepper.error, "slope ", sys()

#Sliding features are on hold until I can write manifold state code.
#sys2 = solver.SlidingDemo()
#sol = solver.Solver(sys2, events=True, slide=True, min_ratio=.2)
#print "\nWith Events, With Sliding - satellite control"
#start = time.clock()
#sys2.x_out, sys2.t_out = sol.simulate(5)
#print "time elapsed ", time.clock() - start
#sys2.plot()

#Meanwhile, brute force works well.
sys3 = lyapunov.SlidingDemo()
sol = lyapunov.Solver(sys3, points=10000)
print "\nNo Events, No Sliding - satellite control"
print "initial state", sys3.state
start = time.clock()
sys3.x_out, sys3.t_out = sol.simulate(5)
print "time elapsed", time.clock() - start
sys3.plot()

sys4 = lyapunov.SubsystemDemo()
sol = lyapunov.Solver(sys4)
print "\nMass-Spring-Damper w/PID control"
print "initial state", sys4.state
start = time.clock()
sys4.x_out, sys4.t_out = sol.simulate(10)
print "time elapsed", time.clock() - start
sys4.plot()


#Need a wrapper object to test lyapunov.Filter (to provide an input signal).
class FilterDemo(object):
	def __init__(self):
		self.plant = lyapunov.Filter((1.0, 3.0, 3.0))
		self.time = 0.0
		self.plant.signal = self.reference
		self.plant.state = (0.0,)*3

	def reference(self):
		return 0.0 if self.time < 1.0 else 1.0
	
	@property
	def state(self):
		return self.plant.state

	@state.setter
	def state(self, x):
		self.plant.state = x

	def __call__(self):
		return self.plant()

	def __len__(self):
		return len(self.plant)

	def plot(self):
		plt.figure()
		try:
			plt.plot(sys5.t_out, sys5.x_out[:,0])
			plt.show()
		except:
			plt.close()


sys5 = FilterDemo()
sol = lyapunov.Solver(sys5)
print "\nFilter"
print "initial state", sys5.state
start = time.clock()
sys5.x_out, sys5.t_out = sol.simulate(10)
print "time elapsed", time.clock() - start
sys5.plot()

