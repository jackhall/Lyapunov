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

import sys
import math
import matplotlib.pyplot as plt
import time
import numpy
from itertools import izip
import pdb
import lyapunov

class Motor(object):
	def __init__(self):
		self.time = 0.0
		self.Rs = 10 	#stator winding resistance - ohms
		self.Ls = 0.14 	#stator winding inductance - henrys
		self.Rr = 10 	#rotor winding resistance - ohms
		self.Lr = 0.06 	#rotor winding inductance - henrys
		self.Vr = 100.4	#rotor voltage - volts
		self.B = 0.01 	#viscous friction - N.m.s/rad
		self.K = 0.714 	#back-emf / torque constant - V.s/rad
		self.J = 0.25 	#rotational inertia - kg.m2
		self.alpha = self.Rs / self.Ls
		self.beta = self.Rr / self.Lr
		self.gamma = self.Vr / self.Lr
		self.a = self.K * self.Ls / self.Lr
		self.b = self.B / self.J
		self.c = self.K * self.Ls / self.J
		self.u = None

	def __len__(self):
		return 4

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, x):
		self._state = x
		self._output = self.h_complete(x)

	def d(self):
		""" Disturbance force at a particular time step.
			Must be function with respect to time. """
		return 0.0 #no disturbance by default

	def output(self):
		return self._output
	
	def __call__(self):
		return self.f(self._state, self.u(), self.d())

	def f(self, x, u, d):
		x1, x2, x3, __ = x
		x1dot = -self.alpha*x1 + u/self.Ls
		x2dot = -self.beta*x2 + self.gamma - self.a*x1*x3
		x3dot = -self.b*x3 + self.c*x1*x2 + d/self.J
		x4dot = x3
		return x1dot, x2dot, x3dot, x4dot

	def h_complete(self, x):
		x1, x2, x3, x4 = x
		y = x4
		ydot = x3
		y2dot = -self.b*x3 + self.c*x1*x2
		return y, ydot, y2dot


class FBLController(object):
	"""feedback linearizing control"""
	def __init__(self, plant):
		self.plant = plant
		self._gains = (1.95, 4.69, 3.75) #(k1, k2, k3)
		self.x = self.r = self.y = None

	def u(self):
		return self._control_effort

	def _u_eq(self):
		epsilon = 0.001
		#Compute the equivalent torque that makes the system linear.
		p = self.plant 
		x1, x2, x3, __ = self.x()
		r, rdot, r2dot, r3dot = self.r()
		#if abs(x2) < epsilon?
		den = epsilon if x2 == 0.0 else p.c*x2
		return (p.Ls * (r3dot + (p.alpha + p.beta)*p.c*x1*x2 - p.gamma*p.c*x1 
					+ p.a*p.c*x3*x1**2 - x3*p.b**2 - p.b*p.c*x1*x2) / den)

	def __call__(self):
		"""takes reference(+derivatives) & states, returns control effort"""
		p = self.plant
		r, rdot, r2dot, r3dot = self.r()
		#Now compute the control torque.
		_, x2, _, _ = self.x()
		den = epsilon if x2 == 0.0 else p.c*x2
		k1, k2, k3 = self._gains
		y, ydot, y2dot = self.y()
		u_c = p.Ls * (k1*(r-y) + k2*(rdot-ydot) + k3*(r2dot-y2dot)) / den
		self._control_effort = self._u_eq() + u_c


class Observer(object):
	def __init__(self, plant):
		self.plant = plant
		self.Lmax = 1000.0 #this value is arbitrary
		self.y = self.u = None

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, xhat):
		self._state = xhat
		x1, x2, x3, x4 = xhat
		#computing observer gains...
		epsilon = 1	#a small number to prevent division by zero
		p = self.plant
		L4 = 25 - p.b - p.beta - p.alpha
		L3 = p.alpha*L4 + p.a*p.b*p.c*p.beta*x1**2 - 234
		den = p.c*(p.a*x1*x3 + x2*(p.alpha+p.beta)) 
		#if abs(den) > epsilon?
		den = den if den != 0.0 else epsilon
		L1 = ((1526 - (L3 + (p.b + p.beta)*L4 + (L4 + (p.alpha - 1))
			 * p.alpha*p.beta*p.a*p.b*p.c*x1**2)*p.alpha**2) / den)
		den = p.c*x1 if x1 != 0.0 else epsilon
		L2 = ((-(p.alpha + p.beta)*L3 + p.alpha*L4*(p.b + p.beta) 
			+ (L4 + p.alpha)*p.beta*p.a*p.b*p.c*x1**2 - p.c*x2*L1 - 997) / den)
		#Set maximum gains! Otherwise they'll blow up when your estimated
		#system loses observability. This is also fixed with epsilon as above.
		L1 = L1 if abs(L1) < self.Lmax else (self.Lmax if L1>0 else -self.Lmax)
		L2 = L2 if abs(L2) < self.Lmax else (self.Lmax if L2>0 else -self.Lmax)
		L3 = L3 if abs(L3) < self.Lmax else (self.Lmax if L3>0 else -self.Lmax)
		self._gains = L1, L2, L3, L4
		#computing estimated output...
		self._output = self.plant.h_complete(xhat)

	def __len__(self):
		return 4

	def output(self):
		return self._output

	def __call__(self):
		xdot = self.plant.f(self._state, self.u(), 0.0)
		output_error = self.y()[0] - self._state[3] 
		return tuple(
			[xidot - output_error*L for (xidot, L) in zip(xdot, self._gains)])


class SMController(FBLController):
	"""sliding mode control"""
	def __init__(self, plant, eta, lmbda):
		self.plant = plant
		self.eta = eta
		self.lmbda = lmbda
		self.x = self.r = self.y = None

	def __call__(self):
		"""needs full reference signal, including derivatives"""
		r, rdot, r2dot, r3dot = self.r()
		#Compute restoring torque
		y, ydot, y2dot = self.y()
		s = (r-y)*self.lmbda**2 + 2*(rdot-ydot)*self.lmbda + (r2dot - y2dot)
		#u_d = self.eta if s > 0 else -self.eta
		u_d = self.eta*math.tanh(10.0*s)
		self._control_effort = self._u_eq() + u_d


#Choose input function
t_in = lyapunov.Time(initial=0, step_size=0.001)
if "step" in sys.argv:
	t_in.final = 8.0
	reference = lyapunov.StepSignal(step_time=2.0, y_final=2.0)
elif "squarewave" in sys.argv:
	t_in.final = 30.0
	reference = lyapunov.SquareWave(period=10.0, y_lower=-1.0, y_upper=1.0)
elif "chirp" in sys.argv:
	t_in.final = 60.0
	reference = lyapunov.ChirpSignal(freq_fcn=lambda t:0.01*t)
else:
	print "defaulting control signal to regulation"
	t_in.final = 8.0
	reference = lyapunov.StepSignal(step_time=t_in.final+1)
t_in._construct()

#Construct subsystems.
plant = Motor()
if "fbl" in sys.argv:
	controller = FBLController(plant)
	print "Feedback Linearizing control"
elif "smc" in sys.argv:
	#does lmbda > 1 make sense?
	controller = SMController(plant, eta=10, lmbda=1.5) 
	print "Sliding Mode control"
else:
	raise RuntimeError("Define a controller! (fbl or smc)")
prefilter = lyapunov.Filter((244, 117.2, 18.75))

#Connect subsystems
plant.u = controller.u
#no disturbance yet		self.plant.d.link_to
controller.r = prefilter.output
prefilter.signal = lambda : reference.value #value is a property
labels = {'reference angle': lambda : reference.value,
		  'filtered angle': lambda : prefilter.state[0],
		  'motor angle': lambda : plant.state[3]}
plant.state = (1.0, 1.0, 0.0, -1.0) #randomize?
prefilter.state = (0.0,)*len(prefilter.state)

if "observe" not in sys.argv:
	controller.y = plant.output
	controller.x = lambda : plant.state  
	system = lyapunov.CompositeSystem([reference, prefilter, 
									controller, plant])
	print "no observer"
else:
	observer = Observer(plant)
	controller.y = observer.output
	controller.x = lambda : observer.state
	observer.y = plant.output
	observer.u = controller.u
	observer.state = (0.3,)*4 #randomize?
	labels['observed angle'] = lambda : observer.state[3]
	system = lyapunov.CompositeSystem([reference, prefilter, 
									controller, plant, observer])
	print "with observer"

#Configure plotter
if "showu" in sys.argv:
	#labels['control effort'] = controller.u
	labels['u_eqivalent'] = controller._u_eq
plotter = lyapunov.Plotter(system, labels)

#Simulate
print "initial state:", system.state
#stepper = solvers.Stepper(sys)
#stepper.step(0.01)
#print "next state:", sys.state
print "simulating for", t_in.span, "sec with", t_in.points, "points."
start = time.clock()
sol = lyapunov.Solver(system, plotter=plotter)
plotter.x, plotter.t = sol.simulate(t_in)
print "final state:", plotter.x[-1]
print "elapsed time =", time.clock() - start

#Plot
plotter.time_response() 
if "showx" in sys.argv:
	print "plotting states..."
	titles = ["stator current", "rotor current", "angular velocity", "angle"]
	for i in range(4):
		plt.figure()
		plt.plot(plotter.t, plotter.x[:,3+i], label="plant")
		if "observe" in sys.argv:
			plt.plot(plotter.t, plotter.x[:,7+i], label="observer")
		plt.legend()
		plt.xlabel('time (s)')
		plt.title(titles[i])
		plt.show()
	print "done"

