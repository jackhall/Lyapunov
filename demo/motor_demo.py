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

import matplotlib.pyplot as plt
import pdb
import lyapunov

###########################
#FBL on motor without observer

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

	def __call__(self):
		"""takes reference(+derivatives) & states, returns control effort"""
		epsilon = 0.001
		#Compute the equivalent torque that makes the system linear.
		p = self.plant 
		x1, x2, x3, __ = self.x()
		r, rdot, r2dot, r3dot = self.r()
		den = epsilon if x2 == 0.0 else p.c*x2
		u_eq = (p.Ls * (r3dot + (p.alpha + p.beta)*p.c*x1*x2 - p.gamma*p.c*x1 
				+ p.a*p.c*x3*x1**2 - x3*p.b**2 - p.b*p.c*x1*x2) / den)
		#Now compute the control torque.
		k1, k2, k3 = self._gains
		y, ydot, y2dot = self.y()
		u_c = p.Ls * (k1*(r-y) + k2*(rdot-ydot) + k3*(r2dot-y2dot)) / den
		self._control_effort = u_eq + u_c


class FBLNoObsrv(object):
	def __init__(self):
		self.time = 0.0
		self.output_labels = ['reference angle', 'filtered reference',
							  'motor angle']
		self.step_time = 0.001
		#construct subsystems
		self.plant = Motor()
		self.controller = FBLController(self.plant)
		self.prefilter = lyapunov.Filter((244, 117.2, 18.75))
		#link statements to connect subsystems
		self.plant.u = self.controller.u
		#no disturbance yet		self.plant.d.link_to
		self.controller.x = lambda : self.plant.state #state is a property
		self.controller.y = self.plant.output
		self.controller.r = self.prefilter.output
		self.prefilter.signal = self.reference
		self.plant.state = (0.0,)*4
		self.prefilter.state = (0.0,)*len(self.prefilter)
		self.controller()

	def __len__(self):
		return len(self.plant) + len(self.prefilter)

	def __call__(self):
		return self.plant() + self.prefilter()

	@property
	def state(self):
		return self.plant.state + self.prefilter.state

	@state.setter
	def state(self, x):
		n = len(self.plant)
		self.plant.state = x[:n]
		self.prefilter.state = x[n:]
		self.controller()

	def reference(self):
		return 0.0 if self.time < self.step_time else 2.0

	@property
	def output(self):
		return (self.reference(), self.prefilter.state[0], self.plant.output()[0])

	def plot(self):
		plt.figure()
		#Use descriptors for labels?
		for i, ilabel in enumerate(self.output_labels):
			plt.plot(self.t_out, self.y_out[:,i], label=ilabel)
		plt.show()

		plt.figure()
		for i in range(4):
			plt.plot(self.t_out, self.x_out[:,i])
		plt.show()

##############################
#FBL with observer

class Observer(object):
	def __init__(self, plant):
		self.plant = plant
		self.state = 0.0, 0.0, 0.0, 0.0
		self.Lmax = 100.0 #this value is arbitrary
		self.y = self.u = None

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, xhat):
		self._state = xhat
		x1, x2, x3, x4 = xhat
		#computing observer gains...
		epsilon = 0.001	#a small number to prevent division by zero
		p = self.plant
		L4 = 25 - p.b - p.beta - p.alpha
		L3 = p.alpha*L4 + p.a*p.b*p.c*p.beta*x1**2 - 234
		L1 = ((1526 - (L3 + (p.b + p.beta)*L4 + (L4 
			+ (p.alpha - 1))*p.alpha*p.beta*p.a*p.b*p.c*x1**2)*p.alpha**2) 
			/ (p.c*(p.a*x1*x3 + x2*(p.alpha+p.beta)) + epsilon))
		L2 = ((-(p.alpha + p.beta)*L3 + p.alpha*L4*(p.b + p.beta) 
			+ (L4 + p.alpha)*p.beta*p.a*p.b*p.c*x1**2 - p.c*x2*L1 - 997) 
			/ (p.c*x1 + epsilon))
		#Set maximum gains! Otherwise they'll blow up when your estimated
		#system loses observability. This is also fixed with epsilon as above.
		#L1 = L1 if L1 < self.Lmax else self.Lmax
		#L2 = L2 if L2 < self.Lmax else self.Lmax
		self._gains = L1, L2, L3, L4
		#computing estimated output...
		self._output = self.plant.h_complete(xhat)

	def output(self):
		return self._output

	def __call__(self):
		xdot = self.plant.f(self._state, self.u(), 0.0)
		output_error = self.y()[0] - self._state[3] 
		return tuple(
			[xidot + output_error*L for (xidot, L) in zip(xdot, self._gains)] )


class FBLObsrv(FBLNoObsrv):
	def __init__(self):
		FBLNoObsrv.__init__(self)
		self.observer = Observer(self.plant)
		self.controller.x = lambda : self.observer.state
		self.observer.y = self.plant.output
		self.observer.u = self.controller.u 
		#disturbance is still assumed to be zero
		self.output_labels.append('estimated angle')

	def __len__(self):
		return len(self.plant) + len(self.prefilter) + len(self.observer)
	
	def __call__(self):
		qdot = self.prefilter()
		xhatdot = self.observer()
		xdot = self.plant() 
		return xdot + qdot + xhatdot

	@FBLNoObsrv.state.getter
	def state(self):
		return self.plant.state + self.prefilter.state + self.observer.state

	@FBLNoObsrv.state.setter
	def state(self, x):
		n = len(self.plant)
		m = len(self.prefilter)
		self.plant.state = x[:n]
		self.prefilter.state = x[n:(n+m)]
		self.observer.state = x[(n+m):]
		self.controller()

	@property
	def output(self):
		return (self.reference(), self.prefilter.state[0], 
				self.plant.output[0], self.observer.output[0]) 

############################
#Sliding mode control - with observer

class SMController(object):
	"""sliding mode control"""
	def __init__(self, plant, eta_value):
		System.__init__(self, plant.state)
		self.eta = eta_value
		self.x = self.r = self.y = None

	@property
	def mode(self):
		"""read self.x and compute mode"""
		pass

	def __call__(self):
		"""needs full reference signal, including derivatives"""
		x1, x2, x3, x4 = self.x()
		u_eq = 0.0
		u_d = 0.0
		return u_eq + u_d


class SMCObsrv(FBLObsrv):
	def __init__(self):
		FBLObsrv.__init__(self)
		self.controller = SMController(self.plant)
		self.plant.u = self.observer.u = self.controller.u
		self.controller.r = self.prefilter.output
		self.controller.y = self.plant.output
		self.controller.x = lambda : self.observer.state

