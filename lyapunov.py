#!/usr/bin/env python
"""
Lyapunov is a toolbox for integrating dynamical systems. 

Instead of treating systems as functions, lyapunov represents systems as 
objects. Not only does this significantly clean up the solver interface, but 
it also encourages the encapsulation of subsystems. Exposed classes:

Recorder - A way to track and later plot any system data. 
CompositeSystem - A container/manager for subsystems, itself a system.
Filter - A linear filter of arbitrary order. 
Time - A convenient way of creating and manipulating time iterables.
StepSignal - Generates a step signal.
SquareWave - Generates a square wave.
SineWave - Generates an sinusoid.
ChirpSignal - Generates a sinusoid with an arbitrary instantaneous frequency.

Exposed functions:

simulate - a way to numerically integrate systems

--For a full description of Lyapunov's system concept, see lyapunov.simulate. 

--lyapunov.CompositeSystem provides a subsystem manager that itself 
  implements the system interface, allowing the user to build arbitrarily 
  complex hierarchies of subsystems.

--For convenient plotting state or output trajectories (in time or as a 
  phase portrait), see lyapunov.Plotter. 

--The file 'demo/motor_demo.py' has a full demonstration of the above 
  features. 

--Event detection is included. See lyapunov.simulate for more information.

--Code generation from symbolic input (from sympy) may happen at some point.

Integration of ordinary equations is done with solvers from 
boost::numeric::odeint, wrapped using boost.python. See the makefile in the 
main directory for tips on compiling. A C++11 capable compiler will be needed, 
along with a recent copy of boost (1.53 or later). Compilation and linking 
will result in a file called 'solvers.so' (or whatever suffix shared 
libraries have on your system). Either place this file in the main directory 
with 'lyapunov.py', or put both in whatever system directory python looks in 
to import external modules.
"""

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
import collections
from itertools import chain, compress, imap, izip, ifilter
import operator
import numpy
import inspect
import matplotlib.pyplot as plt
import solvers

#################
# System manipulation
class ParallelSystems(object):
	"""
	A container/manager for subsystems, itself a system.

	ParallelSystems acts as an aggregate of other objects which implement
	some or all of the system concept. It will ascertain whether each 
	contained object has a 'state' attribute, a 'time' attribute, and/or
	is callable. Objects with state must be callable. ParallelSystems will 
	ignore objects with none of these interface elements. 
	
	The order in which subsystems are stored determines the order of 
	evaluation and update. For instance, when ParallelSystems is called, 
	it will call each of its callable subsystems in order.  When 'state' is 
	queried or set, ParallelSystems will accordingly chain together all subsystem 
	state tuples or slice them.
	"""
	def __init__(self, sys_list):
		"""
		Initialize with an iterable of the subsystems to be managed. 
		Subsystem 'time' attributes are not synchronized until 
		CompositeSystem's 'time' attribute is set.

		Usage: CompositeSystem(sys_list)
		"""
		self._subsystems = list(sys_list) #use an OrderedDict?
		self._dof = [len(sys.state[1]) for sys in sys_list]
		self._time = sys_list[0].state[0] or 0.0
		if not self.are_synchronized():
			print "Subsystem times aren't synchronized yet."
		#Check to make sure that any system that has state is callable.
		for dof, can_call in zip(self._dof, map(callable, sys_list)):
			if dof > 0 and not can_call:
				raise NotImplementedError("Systems with states " 
						+ "should be callable")

	def __call__(self):
		"""
		Call each callable subsystem in turn, and concatenate the results.
		"""
		call_iter = ifilter(callable, self._subsystems)
		return tuple(chain.from_iterable((sys() for sys in call_iter)))

	@property
	def state(self):
		"""
		Concatenate and return state information from all subsystems 
		that have it.
		"""
		return self._time, tuple(chain.from_iterable(
								 (sys.state[1] for sys in self._subsystems)))

	@state.setter
	def state(self, t_x):
		"""
		Distributes slices of a concatenated state tuple to their
		corresponding stated subsystems."""
		self._time, x = t_x
		a = 0
		for dof, sys in izip(self._dof, self._subsystems):
			b = a + dof 
			sys.state = self._time, x[a:b]
			a = b

	def are_synchronized(self, index=None):
		if index is None:
			return not any(imap(lambda sys: sys.state[0] != self._time, 
								self._subsystems))
		else:
			return self._subsystems[index].state[0] == self._time

	def synchronize(self, index=None):
		if index is None:
			for sys in self._subsystems:
				sys.state = self._time, sys.state[1]
		else:
			sys = self._subsystems[index]
			sys.state = self._time, sys.state[1]

	def add_subsystem(self, index, new_system):
		"""
		Starts managing a new subsystem, inserted into the update order
		at index. All subsequent subsystems are shifted back.

		Usage: [CompositeSystem].add_subsystem(index, new_system)
		"""
		self._subsystems.insert(index, new_system)
		self._dof.insert(index, len(new_system.state[1]))
		if not self.are_synchronized(index):
			print "warning - forced to synchronize new subsystem"
			self.synchronize(index)

	def remove_subsystem(self, index):
		"""
		Ceases to manage the subsystem at index in the update order. All
		subsequent subsystems are moved up.

		Usage: [CompositeSubsystem].remove_subsystem(index)
		"""
		self._subsystems.pop(index)
		self._dof.pop(index)


State = collections.namedtuple('State', 't, x')

def state_property(tname='_lyapunov__t', xname='_lyapunov__x'):
	""" 
	Returns a property that acts like a StateTuple. The elements are 
	actually stored in attributes called 'xname' and 'tname' for each
	instance. Usage:

	state(tname='_lyapunov__t', xname='_lyapunov__x') -> property attribute
	"""
	def fget(obj):
		return State(obj.__dict__[tname], obj.__dict__[xname])
	def fset(obj, t_x):
		obj.__dict__[tname], obj.__dict__[xname] = t_x
	def fdel(obj):
		del obj.__dict__[tname], obj.__dict__[xname]
	doc = "Acts like a StateTuple."
	return property(fget, fset, fdel, doc)


def systemfunctor(sys_cls):
	argspec = inspect.getargspec(sys_cls.__call__)
	if len(argspec.args) == len(argspec.defaults)+1:
		#create a new function with two arguments inserted, defaulted to None
		#it may not be possible to do this...
		pass
	else:
		#check varargs?
		raise NotImplementedError("System not callable without arguments")

#################

#################
# Controllers
class PID(object):
	""" SISO """
	def __init__(self, Kp=1.0, Ki=0.0, Kd=0.0):
		self.Kp, self.Ki, self.Kd = Kp, Ki, Kd
		self.r = self.y = None
		self.state = 0.0, (0.0,) #integral term

	state = state_property(xname="_state")

	def __call__(self):
		x, _ = self.y() 
		error = self.r() - x
		return (error,)

	def u(self):
		#catch Nonetype exceptions for y and r!
		x, v = self.y() #assumes second-order, y returns output derivative
		error = self.r() - x
		return self.Kp*error + self.Ki*self.state.x[0] - self.Kd*v

#################

#################
# Input Signals
class StepSignal(object):
	"""
	Generates a step signal. The value of the step signal depends on the
	value of '[StepSignal].time'.
	"""

	def __init__(self, step_time=1.0, y_initial=0.0, y_final=1.0):
		"""
		Initialize with the time at which the step occurs, the value
		before the step, and the value after the step.

		Usage: StepSignal(step_time=1.0, y_initial=0.0, y_final=1.0)
		"""
		self.step_time = step_time
		self.initial, self.final = y_initial, y_final
		self.state = 0.0, ()
	
	state = state_property("time", "_lyapunov__x")

	@property
	def value(self):
		""" The current value of the signal. """
		return self.initial if self.time < self.step_time else self.final


class SquareWave(object):
	"""
	Generates a square wave. The value depends on '[SquareWave].time'.
	"""
	def __init__(self, period=1.0, y_lower=-1.0, y_upper=1.0):
		"""
		Initialize with the period, the trough value, and the peak value.

		Usage: SquareWave(period=1.0, y_lower=-1.0, y_upper=1.0)
		"""
		self.period, self.lower, self.upper = period, y_lower, y_upper
		self.state = 0.0, ()

	state = state_property("time", "_lyapunov__x")

	@property
	def value(self):
		""" The current value of the signal. """
		if self.time % self.period < 0.5*self.period:
			return self.upper
		else:
			return self.lower


class SineWave(object):
	"""
	Generates a sinusoid. Value depends on '[SineWave].time'.
	"""
	def __init__(self, frequency=1.0, mean=0.0, amplitude=1.0, phase=0.0):
		"""
		Initialize with frequency, DC magnitude, AC magnitude, and phase shift.
		Frequency is in rad/s, and the phase shift is positive. 

		Usage: SineWave(frequency=1.0, mean=0.0, amplitude=1.0, phase=0.0)
		"""
		self.frequency, self.phase = frequency, phase
		self.mean, self.amplitude = mean, amplitude
		self.state = 0.0, ()

	state = state_property("time", "_lyapunov__x")
	
	@property
	def value(self):
		""" The current value of the signal. """
		return self.amplitude*math.sin(self.frequency*self.time 
				+ self.phase) + self.mean


class ChirpSignal(object):
	"""
	Generates a sinusoid with an arbitrary instantaneous frequency.
	"""
	def __init__(self, freq_fcn=None, mean=0.0, amplitude=2.0):
		""" 
		Initialize with a function that computes frequency from time, 
		the DC magnitude, the AC magnitude, and a function that computes 
		frequency from time. Frequencies are in rad/s. 

		Usage: ChirpSignal(freq_fcn=None, mean=0.0, amplitude=2.0)
		"""
		self.amplitude, self.mean = amplitude, mean
		self.state = 0.0, ()
		if freq_fcn is None:
			self.freq_fcn = lambda time : time
		else:
			self.freq_fcn = freq_fcn

	state = state_property("time", "_lyapunov__x")

	@property
	def value(self):
		""" The current value of the signal. """
		inst_freq = self.freq_fcn(self.time)
		return self.mean + self.amplitude*math.sin(inst_freq*self.time)

#################

class Filter(object):
	""" 
	Differentiates a reference signal with a linear filter. The 
	'output' attribute refers to the complete reference signal. 
	"""

	def __init__(self, gains):
		""" 
		Place poles before constructing. 'gains' is a list of
		coefficients of the characteristic equation (normalized),
		from lowest-highest order. Exclude the highest, since it 
		should be equal to one anyway. 
		"""
		#This way, there's no need to handle complex numbers.
		self._num_states = len(gains)
		self._gains = gains #check signs?
		self.state = 0.0, (0.0,)*self._num_states #signal isn't connected yet
		self.signal = None

	state = state_property("_time", "_state")

	@state.setter
	def state(self, t_x):
		""" Set state and computes the only nontrivial derivative. """
		self._time, self._state = t_x
		#parenthesis for (q,d)?
		#catch Nonetype exception?
		self._xndot = sum(-q*d for (q,d) in zip(self._state, self._gains)) 

	def output(self):
		return self._state + (self._xndot + self._gains[0]*self.signal(),)

	def __call__(self):
		return self._state[1:] + (self._xndot + self._gains[0]*self.signal(),)


class StopIntegration(Exception):
	""" 
	Exception to be thrown when the system encounters an
	Event and needs to stop integrating.
	"""

	def __iter__(self, name=""):
		self.name = name

	def __str__(self):
		return repr(self.name)


class Time(object):
	"""
	An iterable that acts like a time list.

   	The Time class supports flexible methods for defining said list. Enough 
	of these attributes need to be defined before a given instance
	is used as an iterable (like during a simulation):

	--initial
	--step_size
	--final
	--points (an int)
	--span

	Any of these can be specified as keyword arguments on initialization
	or set as simple data attributes. No checking is performed until the
	'construct' method is called (it's called automatically when the object
	is used as an iterable). See 'construct' docstring for more information.
	"""
	def __init__(self, **kwargs):
		"""
		Accepts any of the following keyword arguments:

		--initial
		--step_size
		--final
		--points (an int)
		--span

		Provided arguments are assigned to the relevant data attributes.
		"""
		self.initial = kwargs.get('initial')
		self.step_size = kwargs.get('step_size')
		self.final = kwargs.get('final')
		self.points = kwargs.get('points')
		self.span = kwargs.get('span')

	def _construct(self):
		"""
		Infer 'initial', 'step_size', and/or 'final'.

		Usage: [Time].construct()

		If any of the primary attributes ('initial', 'final', and
		'step_size') are not provided, the method will attempt to infer them 
		from secondary attributes ('span', 'points'). In case of conflict,
		information from primary attributes is preferred, and secondary 
		attributes will be made consistent.
		"""
		if self.final is None and self.inital_time is None:
			raise RuntimeError("No absolute time information given.")
		if self.step_size is None and self.points is None:
			raise RuntimeError("No time sparseness/density information given.")
		if self.final is None:
			if self.span is not None:
				self.final = self.initial + self.span
			elif self.points is not None and self.dt is not None:
				self.final = (self.initial + 
								   self.step_size * (self.points-1))
			else:
				raise RuntimeError("final undefined")
		if self.initial is None:
			if self.span is not None:
				self.initial = self.final - self.span
			elif self.points is not None and self.dt is not None:
				self.initial = (self.final - 
									 self.step_size * (self.points-1))
			else:
				raise RuntimeError("initial undefined")
		self.span = self.final - self.initial
		if self.step_size is None:
			self.step_size = self.span / self.points
		else:
			self.points = self.span / self.step_size

	def __iter__(self):
		"""
		Use Time objects as an iterable. Note: 'construct' is called.
		"""
		self._construct()
		for t in numpy.linspace(self.initial, self.final, self.points):
			yield t


class Recorder(object):
	"""
	Records system data during a Solver simulation through callbacks.

	Use Plotter to record arbitrary system data during a simulation and
	plot it afterwards. This object is instantiated with a dictionary 
	mapping string variable names to callback functions. These functions
	will be called with no arguments at each iteration of a solving loop 
	to retrieve some scalar value. Using lambdas or system methods is  
	usually convenient.

	After simulation, call either 'time_response' or 'phase_portrait' to
	plot. The former will plot all recorded quantities against time. The
	latter should be called with two labels, which tell it which two 
	variables to plot against each other.

	Calling 'clear' will preserve labels and callbacks, but delete all saved
	variable data.
	"""
	def __init__(self, system, labels={}):
		"""
		Create a Plotter from a dict mapping variable labels to 
		callback functions. The callback functions should be callable
		with no argments and return a scalar numeric. 
		"""
		#labels is a dict ... explain
		self.system = system 
		self.labels = labels 
		self.lines = {label: [] for label in labels.keys()}
		self.x = []
		self.t = []

	def log(self, events=[]):
		"""
		Call with no arguments to record system variables at the 
		current state and time. Usually only used by 'Solver.simulate'.
		"""
		self.t.append(self.system.state[0])
		self.x.append(self.system.state[1])
		for label, f in self.labels.iteritems():
			self.lines[label].append(f())

	def clear(self):
		"""
		Call with no arguments to erase all saved records of state
		variables. The callback functions and the labels themselves 
		remain.
		"""
		self.x = []
		self.t = []
		for label in self.lines:
			self.lines[label] = []

	def time_response(self):
		"""
		Plots all saved records of system variables against time. String
		labels are given in the legend. 
		"""
		plt.figure()
		for label, data in self.lines.iteritems():
			plt.plot(self.t, numpy.array(data), label=label)
		plt.xlabel("time (s)")
		plt.legend()
		plt.show()

	def phase_portrait(self, xlabel, ylabel):
		"""
		Plots the phase space trajectory of two system variables. This
		means that they are plotted against each other for the currently
		saved trajectory. The underlying vector field is NOT plotted, as
		that would involve computing derivatives (which are not saved). 
		Furthermore, the vector field could only be evaluated as a function
		for a second-order autonomous system.
		"""
		#does not evaluate derivatives! explain
		plt.figure()
		plt.plot(numpy.array(self.lines[xlabel]), 
				 numpy.array(self.lines[ylabel]))
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.show()


def simulate(system, time_sequence, **kwargs):
	"""
	An ODE solver function for numerically integrating system objects.

	Usage:
		record = simulate(system, time_sequence, **kwargs)
		record = simulate(system, time_sequence[, logger][, solver]
						  [, events[, tolerance][, event_handler]])

	A time iterable should provide a float for each time step the 
	solver needs to output (not step sizes!) Note that improperly-specified 
	time iterables may cause the solver to diverge from the solution if 
	combined with a fixed-step solver.

	After solving, the system will be reset to its state and time
	before 'simulate' was called.

	An object implements the system concept by providing 'state' (tuple 
	of floats) and 'time' (float) data attributes or properties. Calling 
	a system object with no arguments must return the current state 
	derivatives (tuple of floats). If 'time' is not specified, it will be 
	added and initialized to zero.
	
	State derivatives should obviously correspond one-to-one with states. 
	Setting 'state' and 'time' attributes should uniquely 
	map to state derivatives. This means that calling the system should not 
	visibly mutate the object. 

	The default solver is a Runge-Kutta 4 fixed step solver.
	In the future, variable step solving may be available as well as a wider
	range of solvers. The solver may be specified by passing in the desired
	solver class object with the 'solver' keyword.

	Event detection can be used by passing an iterable of callable objects 
	with the 'events' keyword. Each event in the iterable should take no
	arguments and return a float. An event occurs when the result of one or
	more event functions crosses zero. A rootfinder will then determine the 
	system state and time at the event to within a certain time tolerance. 
	This tolerance can be specified with the 'tolerance' keyword.

	If the events change system state or behavior, or if the set of events
	tracked over the course of simulation will change, you may want to specify
	an event handler with the 'event_handler' keyword. This is a callable
	that takes an iterable of events that have occurred and returns (as an 
	iterable) the next set of events to be tracked. The handler is only called
	when one or more events occur. The events it returns will replace any 
	previous events.
	"""
	#parse arguments and initialize
	if 'solver' in kwargs:
		stepper_class = kwargs['solver']
	else:
		stepper_class = solvers.runge_kutta4
	events_defined = 'events' in kwargs
	if events_defined:
		if 'tolerance' in kwargs:
			stepper = stepper_class(system, time_sequence, kwargs['events'], 
									kwargs['tolerance'])
		else:
			stepper = stepper_class(system, time_sequence, kwargs['events'])
		if 'event_handler' in kwargs:
	 		event_handler = kwargs['event_handler']
		else:
			event_handler = lambda events : events
	else:
		stepper = stepper_class(system, time_sequence)
	if 'logger' in kwargs:
		logger = kwargs['logger']
	else:
		logger = Recorder(system)

	#simulate
	original_state = system.state	
	if events_defined:
		for t, events in stepper:
			#Check to make sure system states have not become invalid.
			if any( map(math.isnan, system.state[1]) ):
				raise ArithmeticError("System state is NaN!")
			logger.log()
			if len(events) > 0:
				stepper.events = event_handler(events)
	else:
		for t in stepper:
			#Check to make sure system states have not become invalid.
			if any( map(math.isnan, system.state[1]) ):
				raise ArithmeticError("System state is NaN!")
			logger.log()
	system.state = original_state
	return logger



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

