#/usr/bin/python

import time
import lyapunov
import solver

sys = solver.SimpleDemo()
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


#sys2 = solver.SlidingDemo()
#sol = solver.Solver(sys2, events=True, slide=True, min_ratio=.2)
#print "\nWith Events, With Sliding - satellite control"
#start = time.clock()
#sys2.x_out, sys2.t_out = sol.simulate(5)
#print "time elapsed ", time.clock() - start
#sys2.plot()

sys3 = solver.SlidingDemo()
sol = solver.Solver(sys3, points=10000)
print "\nWith Events, No Sliding - satellite control"
start = time.clock()
sys3.x_out, sys3.t_out = sol.simulate(5)
print "time elapsed ", time.clock() - start
sys3.plot()

