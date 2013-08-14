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

fil = lyapunov.Filter((1.0, 3.0, 3.0))
ref = lyapunov.StepSignal()
fil.signal = lambda : ref.value
fil.state = (0.0,)*3
sys5 = lyapunov.CompositeSystem([ref, fil])
t_in = numpy.linspace(0.0, 10.0, 100)
print "\nFilter"
print "initial state", sys5.state
start = time.clock()
record = lyapunov.simulate(sys5, t_in)
print "time elapsed", time.clock() - start
x_out, t_out = numpy.array(record.state), numpy.array(record.time)
plt.figure()
try:
	plt.plot(t_out, x_out[:,0])
	plt.show()
except:
	plt.close()

