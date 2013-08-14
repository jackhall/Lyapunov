#!/usr/bin/python

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

