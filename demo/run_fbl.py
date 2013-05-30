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

import pdb
import time
import numpy
import matplotlib.pyplot as plt
from itertools import izip
import lyapunov
#import solvers
import motor_demo


step_time = 0.1
#Construct subsystems.
plant = motor_demo.Motor()
reference = lyapunov.StepSignal(step_time=step_time, yf=2.0)
controller = motor_demo.FBLController(plant)
prefilter = lyapunov.Filter((244, 117.2, 18.75))
#Connect subsystems.
plant.u = controller.u
#no disturbance yet		self.plant.d.link_to
controller.x = lambda : plant.state  #state is a property
controller.y = plant.output
controller.r = prefilter.output
prefilter.signal = lambda : reference.value
plant.state = (1.0,)*4
prefilter.state = (0.0,)*len(prefilter)
sys = lyapunov.CompositeSystem([reference, prefilter, controller, plant])


labels = {'reference angle': lambda : reference.value,
		  'filtered angle': lambda : prefilter.state[0],
		  'motor angle': lambda : plant.state[3]}
plotter = lyapunov.Plotter(sys, labels)


final_time = 8.0
num_points = 1000
print "\nMotor with Feedback Linearization, no observer"
print "initial state:", sys.state
#stepper = solvers.Stepper(sys)
#stepper.step(0.01)
#print "next state:", sys.state
print "simulating for", final_time, "sec with", num_points, "points."
start = time.clock()
sol = lyapunov.Solver(sys, points=num_points, plotter=plotter)
plotter.x, plotter.t = sol.simulate(final_time)
print "final state:", sol.x_out[-1]
print "elapsed time =", time.clock() - start
plotter.time_response()

plt.figure()
plt.plot(numpy.array(plotter.t), numpy.array(plotter.x)[:,-1])
plt.show()

pdb.set_trace()
