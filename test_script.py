#!/usr/bin/python

import time
import lyapunov
import solver

sys = solver.DemoNoEvents()
sys.time = 0.0
stepper = lyapunov.Stepper(sys)
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


sys2 = solver.DemoEvents()
sol = solver.Solver(sys2, events=True, slide=True, min_ratio=.5)
print "\nWith Events, With Sliding - satellite control"
start = time.clock()
sys2.x_out, sys2.t_out = sol.simulate(1.7)
print "time elapsed ", time.clock() - start
sys2.plot()

#sys3 = solver.DemoEvents()
#sol = solver.Solver(sys3, events=True, slide=False, min_ratio=.1)
#print "\nWith Events, No Sliding - satellite control"
#start = time.clock()
#sys3.x_out, sys3.t_out = sol.simulate(2)
#print "time elapsed ", time.clock() - start
#sys3.plot()

