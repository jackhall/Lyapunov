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
#import solvers
import motor_demo

final_time = 8.0
num_points = 1000
sys = motor_demo.FBLNoObsrv()
sys.state = (1.0,)*4+ (0.0,)*3 #motor + prefilter
sys.time = 0.0
sys.step_time = 0.1
print "\nMotor with Feedback Linearization, no observer"
print "initial state:", sys.state

#stepper = solvers.Stepper(sys)
#stepper.step(0.01)
#print "next state:", sys.state

print "simulating for", final_time, "sec with", num_points, "points."
start = time.clock()
sol = lyapunov.Solver(sys, points=num_points)
sys.x_out, sys.y_out, sys.t_out = sol.simulate(final_time)
print "final state:", sol.x_out[-1]
print "elapsed time =", time.clock() - start
sys.plot()

