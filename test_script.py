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
import numpy
import matplotlib.pyplot as plt

sys = lyapunov.SimpleDemo()
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

#Brute force works well.
#sys3 = lyapunov.SlidingDemo()
##sol = lyapunov.Solver(sys3)
#t_in = numpy.linspace(0.0001, 3.0, 10000)
#record = lyapunov.Recorder(sys3)
#stepper = solvers.runge_kutta4(sys3, t_in)
#print "\nNo Events, No Sliding - satellite control"
#print "initial state", sys3.state
#start = time.clock()
##record = sol.simulate(t_in)
#for t in stepper:
#	record.record()
#print "time elapsed", time.clock() - start
#sys3.x_out, sys3.t_out = numpy.array(record.state), numpy.array(record.time)
#sys3.plot()

#With events and sliding is faster.
sys2 = lyapunov.SlidingDemo()
sys2.u = lambda: -1 #replaces u_naive
t_in = numpy.linspace(0.1, 3.0, 30)
record = lyapunov.Recorder(sys2)
stepper = solvers.runge_kutta4(sys2, t_in, [sys2.s], 0.0001)
print "\nWith Events and Sliding - satellite control"
start = time.clock()
for t, events in stepper:
	record.record()
	if len(events) > 0:
		print "event! ", events
		if events[0] is sys2.s and sys2.u_margin() > 0:
			sys2.u = sys2.u_effective
			stepper.events = [sys2.u_margin]
		else:
			stepper.step_through()
			sys2.u = (lambda: -1) if sys.s() > 0 else (lambda: 1)
			stepper.events = [sys2.s]
print "time elapsed ", time.clock() - start
sys2.x_out, sys2.t_out = numpy.array(record.state), numpy.array(record.time)
sys2.plot()

#sys4 = lyapunov.SubsystemDemo()
#sol = lyapunov.Solver(sys4)
#t_in.step_size, t_in.points = None, 200
#print "\nMass-Spring-Damper w/PID control"
#print "initial state", sys4.state
#start = time.clock()
#record = sol.simulate(t_in)
#sys4.x_out, sys4.t_out = numpy.array(record.state), numpy.array(record.time)
#print "time elapsed", time.clock() - start
#sys4.plot()
#
#
#fil = lyapunov.Filter((1.0, 3.0, 3.0))
#ref = lyapunov.StepSignal()
#fil.signal = lambda : ref.value
#fil.state = (0.0,)*3
#sys5 = lyapunov.CompositeSystem([ref, fil])
#sol = lyapunov.Solver(sys5)
#t_in.step_size, t_in.final = None, 10
#print "\nFilter"
#print "initial state", sys5.state
#start = time.clock()
#record = sol.simulate(t_in)
#x_out, t_out = numpy.array(record.state), numpy.array(record.time)
#print "time elapsed", time.clock() - start
#plt.figure()
#try:
#	plt.plot(t_out, x_out[:,0])
#	plt.show()
#except:
#	plt.close()

