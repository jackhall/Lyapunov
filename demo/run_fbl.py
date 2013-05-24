#!/usr/bin/python
import time
import lyapunov
import motor_demo

final_time = 0.003
num_points = 1000
sys = motor_demo.FBLNoObsrv()
sys.state = (1.0,)*4+ (0.0,)*3 #motor + prefilter
sys.time = 0.0
print "Initial State:", sys.state
print "Simulating for", final_time, "sec with", num_points, "points."
start = time.clock()
sol = lyapunov.Solver(sys, points=num_points)
print "Elapsed time =", time.clock() - start
sys.x_out, sys.y_out, sys.t_out = sol.simulate(final_time)
sys.plot()

