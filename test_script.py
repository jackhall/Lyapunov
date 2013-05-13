#!/usr/bin/python

import lyapunov
import solver

sys = solver.DemoNoEvents()
sys.time = 0.0
stepper = lyapunov.Stepper(sys)
stepper.step(0.1)
print sys.state, sys.time

