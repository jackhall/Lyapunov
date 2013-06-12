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

import math
import pdb
from itertools import chain, compress, imap
import operator
import numpy
import matplotlib.pyplot as plt
import solvers

class SimpleDemo(object):
	"""Mass spring damper system."""
	def __init__(self):
		self.state = (1.0, 1.0) 
		self.u = None

	def __len__(self):
		return 2

	def __call__(self):
		_, v = self._state
		return (v, self.a)

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, x):
		self._state = x
		try:
			self.a = -x[1] - x[0] + self.u()
		except (AttributeError, TypeError): 
			self.a = -x[1] - x[0]

	def plot(self):
		plt.figure()
		try:
			plt.plot(self.x_out[:,0], self.x_out[:,1])
			plt.show()
		except AttributeError:
			plt.close()


class SlidingDemo(object):
	"""Double integrator with linear switching mode."""
	def __init__(self):
		self.state = (1.0, 1.0)
	
	@property
	def mode(self):
		x, v = self._state
		return 'under' if v < -0.5*x else 'over'

	def __len__(self):
		return 2

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, x):
		self._state = x
		self.u = 1 if self.mode == 'under' else -1

	def __call__(self):
		_, v = self.state
		return (v, self.u)

	def plot(self):
		plt.figure()
		try:
			xmin, xmax = min(self.x_out[:,0]), max(self.x_out[:,0])
			ymin, ymax = min(self.x_out[:,1]), max(self.x_out[:,1])
			x = numpy.array([xmin, xmax])
			lam = numpy.array([-0.5*xmin, -0.5*xmax])
			plt.plot(x, lam)
			plt.plot(self.x_out[:,0], self.x_out[:,1])
			plt.show()
		except AttributeError:
			plt.close()


#################
# Controllers
class PID(object):
	def __init__(self, Kp=1.0, Ki=0.0, Kd=0.0):
		self.Kp, self.Ki, self.Kd = Kp, Ki, Kd
		self.r = self.y = None
		self.state = (0.0,) #integral term
	
	def __call__(self):
		x, v = self.y()
		error = self.r() - x
		self._force = self.Kp*error + self.Ki*self.state[0] - self.Kd*v
		return (error,)

	def u(self):
		return self._force

################

class CompositeSystem(object):
	def __init__(self, sys_list):
		self._subsystems = list(sys_list) #just use an OrderedDict!
		self._have_state = map(self._has_state, sys_list)
		self._have_time = map(self._has_time, sys_list)
		self._are_callable = [hasattr(sys, "__call__") for sys in sys_list]
		self._num_states = sum(len(sys) for sys in 
							compress(sys_list, self._have_state))
		self.time = 0.0
		#Check to make sure that any system that has state is callable.
		for has_state, can_call in zip(self._have_state, self._are_callable):
			if has_state and not can_call:
				raise NotImplementedError("Systems with states " 
						+ "should be callable")

	def __call__(self):
		call_iter = compress(self._subsystems, self._are_callable)
		#NoneType error when call_iter is used?
		return tuple(chain.from_iterable(
					 imap(self._eval_systems, call_iter)))

	@staticmethod
	def _eval_systems(sys):
		"""For calling systems."""
		derivative = sys()
		return derivative if derivative is not None else ()

	@property
	def state(self):
		state_iter = compress(self._subsystems, self._have_state)
		return tuple(chain.from_iterable(
					 imap(lambda sys : sys.state, state_iter)))

	@property
	def subsystems(self):
		return list(self._subsystems)

	@state.setter
	def state(self, x):
		a = 0
		for sys in compress(self._subsystems, self._have_state):
			b = a + len(sys)
			sys.state = x[a:b]
			a = b
		assert b == self._num_states

	@property
	def time(self):
		return self._time

	@time.setter
	def time(self, t):
		self._time = t
		for sys in compress(self._subsystems, self._have_time):
			sys.time = t

	def __len__(self):
		return self._num_states
		
	@staticmethod
	def _has_state(sys):
		try:
			sys.state
			return True
		except AttributeError:
			return False	

	@staticmethod
	def _has_time(sys):
		try:
			sys.time
			return True
		except AttributeError:
			return False

	def replace_subsystem(self, index, new_system):
		self._subsystems[index] = new_system
		self._have_state[index] = self._has_state(new_system)
		self._have_time[index] = self._has_time(new_system)
		self._are_callable[index] = hasattr(new_system, "__call__")

	def add_subsystem(self, index, new_system):
		self._subsystems.insert(index, new_system)
		self._have_state.insert(index, self._has_state(new_system))
		self._have_time.insert(index, self._has_time(new_system))
		self._are_callable.insert(index, hasattr(new_system, "__call__"))

	def remove_subsystem(self, index):
		self._subsystems.pop(index)
		self._have_state.pop(index)
		self._have_time.pop(index)
		self._are_callable.pop(index)


class SubsystemDemo(object):
	def __init__(self):
		self.plant = SimpleDemo()
		self.control = PID(Ki=1)
		self.control.y = lambda : self.plant.state
		self.control.r = self.reference
		self.plant.u = self.control.u
		self.time = 0.0

	@property
	def state(self):
		return self.plant.state + self.control.state

	@state.setter
	def state(self, x):
		self.plant.state = x[:-1]
		self.control.state = (x[-1],)

	def __call__(self):
		return self.plant() + self.control()

	def reference(self):
		return 0.0 if self.time < 2.0 else 1.0

	def __len__(self):
		return 3

	def plot(self):
		self.plant.t_out = self.t_out
		self.plant.x_out = self.x_out
		self.plant.plot()
		

#################
# Input Signals
class StepSignal(object):
	def __init__(self, step_time=1.0, y0=0.0, yf=1.0):
		self.step_time, self.initial, self.final = step_time, y0, yf
		self.time = 0.0

	@property
	def value(self):
		return self.initial if self.time < self.step_time else self.final


class SquareWave(object):
	def __init__(self, period=1.0, y_lower=-1.0, y_upper=1.0):
		"""frequency is measured in Hz"""
		self.period, self.lower, self.upper = period, y_lower, y_upper
		self.time = 0.0

	@property
	def value(self):
		if self.time % self.period < 0.5*self.period:
			return self.upper
		else:
			return self.lower


class SineWave(object):
	def __init__(self, frequency=1.0, mean=0.0, amplitude=1.0, phase=0.0):
		""" Frequency is in rad/s. """
		self.frequency, self.phase = frequency, phase
		self.mean, self.amplitude = mean, amplitude
		self.time = 0.0
	
	@property
	def value(self):
		return self.amplitude*math.sin(self.frequency*self.time 
				+ self.phase) + self.mean


class ChirpSignal(object):
	def __init__(self, f0, freq_fcn=None, amplitude=2.0, mean=0.0):
		""" Frequencies are in rad/s. Instantaneous frequency is 
			f = f0 + freq_fcn(time) """
		self.amplitude, self.mean = amplitude, mean
		self.f0 = f0
		self.time = 0.0
		if freq_fcn is None:
			self.freq_fcn = lambda time : time
		else:
			self.freq_fcn = freq_fcn

	@property
	def value(self):
		inst_freq = self.f0 * self.freq_fcn(self.time)
		return self.mean + self.amplitude*math.sin(inst_freq*self.time)

#################

class Filter(object):
	""" Differentiates a reference signal with a linear filter. The 
		'output' attribute refers to the complete reference signal. """
	def __init__(self, gains):
		""" Place poles before constructing. 'gains' is a list of
			coefficients of the characteristic equation (normalized),
			from lowest-highest order. Exclude the highest, since it 
			should be equal to one anyway. """
		#This way, there's no need to handle complex numbers.
		self._num_states = len(gains)
		self._gains = gains #check signs?
		self._state = (0.0,)*(self._num_states) #signal isn't connected yet
		self.signal = None

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, x):
		self._state = x
		#Compute the only nontrivial deriviative.
		self._xndot = ( sum(-q*d for (q,d) in zip(self._state, self._gains)) 
					  + self._gains[0]*self.signal()) 

	def __len__(self):
		return self._num_states

	def output(self):
		return self._state + (self._xndot,)

	def __call__(self):
		return self._state[1:] + (self._xndot,)


class Event(object):
	""" A terminal event - stops solving before a sign change """
	def __init__(self, name, function):
		self.name = name 
		self.function = function

	def __call__(self):
		return self.function()

	@property
	def flag(self):
		return False

	@flag.setter
	def flag(self, x):
		if x is True:
			raise RuntimeError(self.name)


class Mode(object):
	""" A persistent mode, with entry and exit functions """
	def __init__(self, name, enter, exit, active=False):
		""" enter and exit should take no arguments and each return a float
		    that changes sign when an event occurs	"""
		self.name = name
		self.enter, self.exit = enter, exit
		self.flag = active 
		#self.terminal = terminal #whether or not solver halts

	def __call__(self):
		if self.flag:
			return self.exit()
		else:
			return self.enter()

#Event and Mode represent terminal and sticking events, respectively.
#Remember to flag an event as active just before stepping through the boundary!
#Make sure to catch the error properly in order to return state history!
#A common interface for setting system mode is needed. Must work with 
#	CompositeSystem and be reasonably easy to emulate.

class Solver(object):
	def __init__(self, system, events=[], min_ratio=.01, 
			     points=100, plotter=None): 
		self.system = system
		self.plotter = plotter #has a update()
		self.events = events #have enter() and exit(), which return floats
		self.stepper = solvers.Stepper(system)
		#Check basic requirements of a system object...
		system.state #to raise an exception if state is not an attribute
		if not hasattr(system, '__call__'):
			raise AttributeError("Need to compute state derivatives")
		self.points = points #number of points recorded/plotted
		self.min_ratio = min_ratio

	@property
	def min_ratio(self):
		return self._min_step_ratio

	@min_ratio.setter
	def min_ratio(self, new_min_ratio):
		if new_min_ratio > 1:
			self._min_step_ratio = 1.0 / min_ratio
		else:
			self._min_step_ratio = new_min_ratio

	def simulate(self, final_time):
		#Input checks for events an system.time
		events_provided = len(self.events) != 0
		if events_provided:
			active = []
			inactive = self.events
		autonomous = not hasattr(system, 'time')
		if autonomous:
			self.system.time = 0.0 #means that system.time is reserved!
		#Initialize solving loop.
		step_size = (final_time - self.system.time) / self.points
		x_out = [self.system.state]
		t_out = [self.system.time]
		#main solver loop
		while self.system.time < final_time:
			#Step forward in time.
			self.stepper.step(step_size)
			#Check to make sure system states have not become invalid.
			if True in map(math.isnan, self.system.state):
				print "System state is NaN!"
				pdb.set_trace()
			#Detect an event - defined as a change in sign.
			if events_provided:
				if self.system.mode != current_mode:
					self.stepper.find_root(step_size, 
										   step_size*self.min_ratio);
					current_mode = str(self.system.mode)
			#Record for output.
			if self.plotter is not None:
				self.plotter.update()
			x_out.append(self.system.state)
			t_out.append(self.system.time)
		#Reset system to initial conditions.
		self.system.state, self.system.time = x_out[0], t_out[0]
		if autonomous:
		 	del self.system.time
		return numpy.array(x_out), numpy.array(t_out)
			

class Plotter(object):
	def __init__(self, system, labels):
		self.system = system
		self.labels = labels
		self.lines = {label: [] for label in labels.iterkeys()}
		self.time = []

	def update(self):
		self.time.append(self.system.time)
		for label, f in self.labels.iteritems():
			self.lines[label].append(f())

	def clear(self):
		self.time = []
		for label in self.lines:
			self.lines[label] = []

	def time_response(self):
		plt.figure()
		for label, data in self.lines.iteritems():
			plt.plot(self.time, numpy.array(data), label=label)
		plt.xlabel("time (s)")
		plt.legend()
		plt.show()

	def phase_portrait(self, xlabel, ylabel):
		plt.figure()
		plt.plot(numpy.array(self.lines[xlabel]), 
				 numpy.array(self.lines[ylabel]))
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.show()

	def reconstruct(self, solver):
		self.clear()
		for x, t in zip(solver.x_out, solver.t_out):
			self.system.state = x
			self.system.time = t
			self.update()		


class PlotterList(object):
	def __init__(self, plt_list):
		self.plotters = plt_list

	def update(self):
		for plotter in self.plotters:
			plotter.update()


	#def phase_portrait(self):
	#	#Number of arrows per dimension=20 is arbitrary but works well.
	#	assert len(self.system) == 2
	#	x, y = numpy.meshgrid(numpy.linspace(self.xmin, self.xmax, 20), 
	#						  numpy.linspace(self.ymin, self.ymax, 20)) 
	#	Dx, Dy = numpy.array(x), numpy.array(y)
	#	for i in range(x.shape[0]):
	#		for j in range(x.shape[1]):
	#			coordinates = [float(x[i,j]), float(y[i,j])]
	#			#only works for 2D systems so far...
	#			self.system.state = coordinates #is it ok to assume a time slice is ok?
	#			Dx[i,j], Dy[i,j] = self.system()
	#	plt.quiver(x, y, Dx, Dy)

#Here is a pythonic enum pattern:
#def enum(**enums):
#	return type('Enum', (), enums)
#Numbers = enum(ONE=1, TWO=2.2222, THREE="three")
#Numbers.ONE
#Numbers.TWO
#Numbers.THREE

#...and another one
#DIRECTIONS = set(['up', 'down', 'left', 'right'])
#def move(self, direction):
#	# only if you feel like checking
#	assert direction in DIRECTIONS
#	# you can still just use the strings!
#	if direction == 'up':
#		# Do something

