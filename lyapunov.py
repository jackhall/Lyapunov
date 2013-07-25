#!/usr/bin/env python
"""
Lyapunov is a toolbox for integrating dynamical systems. 

Instead of treating systems as functions, lyapunov represents systems as 
objects. Not only does this significantly clean up the solver interface, but 
it also encourages the encapsulation of subsystems. Exposed classes:

Solver - An object-oriented ODE solver.
Plotter - A way to track and later plot any system data. 
CompositeSystem - A container/manager for subsystems, itself a system.
Filter - A linear filter of arbitrary order. 
StepSignal - Generates a step signal.
SquareWave - Generates a square wave.
SineWave - Generates an sinusoid.
ChirpSignal - Generates a sinusoid with an arbitrary instantaneous frequency.

--For a full description of Lyapunov's system concept, see lyapunov.Solver. 

--lyapunov.CompositeSystem provides a subsystem manager that itself 
  implements the system interface, allowing the user to build arbitrarily 
  complex hierarchies of subsystems.

--For convenient plotting state or output trajectories (in time or as a 
  phase portrait), see lyapunov.Plotter. 

--The file 'demo/motor_demo.py' has a full demonstration of the above 
  features. 

--Event detection and constraint management are in the works.

--Code generation from symbolic input (from sympy) may happen at some point.

Integration of ordinary equations is done with a custom solver implemented
in C++ using boost.python. See the makefile in the main directory for tips
on compiling. A C++11 capable compiler will be needed, along with a recent
copy of boost. Compilation and linking will result in a file called 
'solvers.so' (or whatever suffix shared libraries have on your system). 
Either place this file in the main directory with 'lyapunov.py', or put both
in whatever system directory python looks in to import external modules.
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

	#plot functions are deprecated!
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
	""" SISO """
	#need __len__
	def __init__(self, Kp=1.0, Ki=0.0, Kd=0.0):
		self.Kp, self.Ki, self.Kd = Kp, Ki, Kd
		self.r = self.y = None
		self.state = (0.0,) #integral term

	def __call__(self):
		#catch Nonetype exceptions for y and r!
		x, v = self.y() #assumes second-order, y returns output derivative
		error = self.r() - x
		self._force = self.Kp*error + self.Ki*self.state[0] - self.Kd*v
		return (error,)

	def u(self):
		return self._force

################

class CompositeSystem(object):
	"""
	A container/manager for subsystems, itself a system.

	CompositeSystem acts as an aggregate of other objects which implement
	some or all of the system concept. It will ascertain whether each 
	contained object has a 'state' attribute, a 'time' attribute, and/or
	is callable. Objects with state must be callable. CompositeSystem will 
	ignore objects with none of these interface elements. 
	
	The order in which subsystems are stored determines the order of 
	evaluation and update. For instance, when CompositeSystem is called, 
	it will call each of its callable subsystems in order.  When 'state' is 
	queried or set, CompositeSystem will accordingly chain together all subsystem 
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
		self._dof = map(self._count_states, sys_list)
		self._have_state = [n>0 for n in self._dof]
		self._have_time = map(self._has_time, sys_list)
		self._are_callable = [hasattr(sys, "__call__") for sys in sys_list]
		self._num_states = sum(len(sys.state) for sys in 
							compress(sys_list, self._have_state))
		if sum(self._have_time):
			sys_iter = compress(sys_list, self._have_time)
			self._time = sys_iter.next().time
			for sys in sys_iter:
				if sys.time != self._time:
					print "Subsystem times aren't synchronized yet."
					break
		else:
			self._time = 0
		#Check to make sure that any system that has state is callable.
		for has_state, can_call in zip(self._have_state, self._are_callable):
			if has_state and not can_call:
				raise NotImplementedError("Systems with states " 
						+ "should be callable")

	def __call__(self):
		"""
		Call each callable subsystem in turn, and concatenate the results.
		"""
		#generator that skips non-callable subsystems
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
		"""
		Concatenate and return state information from all subsystems 
		that have it.
		"""
		#generator skipping non-state subsystems
		state_iter = compress(self._subsystems, self._have_state)
		return tuple(chain.from_iterable(
					 imap(lambda sys : sys.state, state_iter)))

	@state.setter
	def state(self, x):
		"""
		Distributes slices of a concatenated state tuple to their
		corresponding stated subsystems."""
		a = 0
		state_iter = compress(self._subsystems, self._have_state) 
		dof_iter = compress(self._dof, self._have_state)
		for dof, sys in zip(dof_iter, state_iter):
			b = a + dof 
			sys.state = x[a:b]
			a = b
		assert b == self._num_states

	@property
	def time(self):
		return self._time

	@time.setter
	def time(self, t):
		"""Updates the 'time' attribute for all subsystems that have it."""
		self._time = t
		for sys in compress(self._subsystems, self._have_time):
			sys.time = t

	@staticmethod
	def _count_states(sys):
		try:
			return len(sys.state)
		except AttributeError:
			return 0	

	@staticmethod
	def _has_time(sys):
		try:
			sys.time
			return True
		except AttributeError:
			return False

	def add_subsystem(self, index, new_system):
		"""
		Starts managing a new subsystem, inserted into the update order
		at index. All subsequent subsystems are shifted back.

		Usage: [CompositeSystem].add_subsystem(index, new_system)
		"""
		self._subsystems.insert(index, new_system)
		self._dof.insert(index, self._count_states(new_system))
		self._have_state.insert(index, self._dof[index]>0)
		self._have_time.insert(index, self._has_time(new_system))
		self._are_callable.insert(index, hasattr(new_system, "__call__"))

	def remove_subsystem(self, index):
		"""
		Ceases to manage the subsystem at index in the update order. All
		subsequent subsystems are moved up.

		Usage: [CompositeSubsystem].remove_subsystem(index)
		"""
		self._subsystems.pop(index)
		self._dof.pop(index)
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
		self.time = 0.0

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
		self.time = 0.0

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
		self.time = 0.0
	
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
		self.time = 0.0
		if freq_fcn is None:
			self.freq_fcn = lambda time : time
		else:
			self.freq_fcn = freq_fcn

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
		self._state = (0.0,)*(self._num_states) #signal isn't connected yet
		self.signal = None

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, x):
		""" Set state and computes the only nontrivial derivative. """
		self._state = x
		#parenthesis for (q,d)?
		#catch Nonetype exception?
		self._xndot = ( sum(-q*d for (q,d) in zip(self._state, self._gains)) 
					  + self._gains[0]*self.signal()) 

	def output(self):
		return self._state + (self._xndot,)

	def __call__(self):
		return self._state[1:] + (self._xndot,)


class StopIntegration(Exception):
	""" 
	Exception to be thrown when the system encounters an
	Event and needs to stop integrating.
	"""

	def __iter__(self, name=""):
		self.name = name

	def __str__(self):
		return repr(self.name)


#require __call__() for event function
#require flag() method to return new event (or None) or raise StopIntegration
#associate events with a given system? no need with chaining
#Remember to flag an event as active just before stepping through the boundary!
#Make sure to catch the error properly in order to return state history!

class Event(type):
	""" A terminal event - stops solving before a sign change """
	def __init__(self, function, name="event"):
		self.name = name 
		self.__call__ = function

	@property
	def flag(self):
		return False

	@flag.setter
	def flag(self, x):
		if x is True:
			raise StopSolving(self.name)

#don't mix function and class decorators!
#with above design, decorators may not be needed
def event(function):
	class _event(object):
		def __init__(self, system):
			self.system = system
		def __call__(self):
			return function(self.system)
	return _event(function, function.__name__)


class Mode(object):
	""" A persistent mode, with entry and exit functions """
	def __init__(self, entry, **kwargs):
		""" enter and exit should take no arguments and each return a float
		    that changes sign when an event occurs	"""
		self.fentry = entry
		nil = lambda: None
		self.fapply = kwargs['apply'] if 'apply' in kwargs else nil
		self.fexit = kwargs['exit'] if 'exit' in kwargs else entry
		self.fignore = kwargs['ignore'] if 'ignore' in kwargs else nil
		self._flag = kwargs['active'] if 'active' in kwargs else False
		self.__call__ = exit if self._flag else entry

	@property
	def flag(self):
		return self._flag

	@flag.setter
	def flag(self, x):
		if self._flag:
			assert not x
			self.fignore()
			self.__call__ = self.fentry
		else:
			assert x
			self.fapply()
			self.__call__ = self.fexit
		self._flag = x

	def apply(function):
		self.fapply = function

	def exit(function):
		self.fexit = function

	def ignore(function):
		self.fignore = function

def mode(function):
	return Mode(function)
	

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


class Solver(object):
	"""
	An ODE solver object for numerically integrating system objects.

	Solver integrates the system with which it is initialized. An object
	implements the system concept by providing 'state' (tuple of floats) 
	and 'time' (float) data attributes or properties. Calling a system 
	object with no arguments must return the current state derivatives 
	(tuple of floats). 
	
	State derivatives should obviously correspond one-to-one with states. 
	Setting 'state' and 'time' attributes should completely and uniquely 
	map to state derivatives. This means that calling the system should not 
	change the object in any way visible to Lyapunov). The 'time' attribute 
	is optional (not needed for autonomous systems), but that name is still 
	reserved. 

	Initialize Solver with a system and (optionally) a plotter object. To
	run a simulation, call 'simulate'.
	Information on step size and simulation length should be provided when
	'simulate' is called or beforehand. The state and time of the system
	when 'simulate' is called will be used as the initial state and time.

	The solver is a Cash-Karp Runge-Kutta solver (order 5) with a fixed step.
	In the future, variable step solving may be available as well as a wider
	range of solvers. 

	Event detection is not finished at the moment, but it's on the way.
	Should events be associated with system or solver?
	"""

	def __init__(self, system, plotter=None, events=[]):
		"""
		Instantiate an ODE solver for a given system.

		Usage: Solver(system, time=Time(), plotter=None)

		A minimal system object is required to instantiate a Solver. 
		AttributeErrors will be raised if the system concept is incomplete.

		A plotter object (anything that provides an 'update' method that 
		takes no arguments) is optional for simulation, but can also be 
		provided here.
		"""
		#events and rootfinding options omitted for now
		#is there a way to make plotter modular wrt system?
		self.system = system #has state and __call__() [time optional]
		self.plotter = plotter #has a update()
		self.events = events #each has __call__() and flag()
		self.stepper = solvers.Stepper(system) #where does Stepper come from?
		#Check basic requirements of a system object...
		system.state #to raise an exception if state is not an attribute
		if not hasattr(system, '__call__'):
			raise AttributeError("Need to compute state derivatives")


	def simulate(self, time, initial_state=None):
		"""
		Numerically integrates the system over the given time iterable.

		Such an iterable should provide a float for each time step the 
		solver needs to output. Note that improperly-specified 'time' 
		iterables may cause the solver to diverge from the solution if 
		combined with a fixed-step solver.

		After solving, the system will be reset to its state and time
		before 'simulate' was called.
		"""
		#Initialize system
		try:
			self.system.time
		except AttributeError:
			autonomous = True
		else:
			autonomous = False
			original_time = self.system.time
		#need code here that works with generators AND containers
		time = iter(time)
		self.system.time = time.next() #means that system.time is reserved!
		original_state = self.system.state
		if initial_state is not None:
			self.system.state = initial_state
		#Configure for events
		events_provided = len(self.events) > 0
		sign_change = lambda f, e: f()*e < 0 #only used if events provided
		#Initialize recording
		class Recorder(object):
			def __init__(self, system, plotter):
				self.x_out = [system.state]
				self.t_out = [system.time]
				self.plotter = plotter

			def __call__(self, system):
				self.x_out.append(system.state)
				self.t_out.append(system.time)
				if self.plotter is not None:
					self.plotter.update()

		recorder = Recorder()
		recorder(self.system)

		#main solver loop (might be cleaner as a recursion)
		for t in time:
			#Step forward in time.
			self.stepper.step(t - self.system.time) 
			#Check to make sure system states have not become invalid.
			if True in map(math.isnan, self.system.state):
				print "System state is NaN!"
				pdb.set_trace()
			#Detect an event - defined as a change in sign.
			if events_provided:
				occurred = imap(sign_change, self.events, event_values)
				active_events = compress(self.events, occurred)
				if len(active_events) > 0:
					try:
						#find_root should update self.events
						self.stepper.find_root(self.events)
					except StopIntegration:
						#clean up after simulation and exit
						recorder(self.system)
						break
				event_values = [f() for f in self.events] #for next step
			#Record state and output information
			recorder(self.system)
		#Reset system to original conditions (before 'simulate' was called)
		self.system.state = original_state
		if autonomous:
		 	del self.system.time
		else:
			self.system.time = original_time
		return numpy.array(recorder.x_out), numpy.array(recorder.t_out)
			

class Plotter(object):
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
	#Update to use object-oriented interface from matplotlib. 
	#Use 3-tiered dict to store labels: figures, subfigures, lines?
	#Flat is better than nested. 
	def __init__(self, system, labels):
		"""
		Create a Plotter from a dict mapping variable labels to 
		callback functions. The callback functions should be callable
		with no argments and return a scalar numeric. 
		"""
		#labels is a dict ... explain
		self.system = system #only used for time & reconstruct
		self.labels = labels 
		self.lines = {label: [] for label in labels.iterkeys()}
		self.time = []

	def update(self):
		"""
		Call with no arguments to record system variables at the 
		current state and time. Usually only used by 'Solver.simulate'.
		"""
		self.time.append(self.system.time)
		for label, f in self.labels.iteritems():
			self.lines[label].append(f())

	def clear(self):
		"""
		Call with no arguments to erase all saved records of state
		variables. The callback functions and the labels themselves 
		remain.
		"""
		self.time = []
		for label in self.lines:
			self.lines[label] = []

	def time_response(self):
		"""
		Plots all saved records of system variables against time. String
		labels are given in the legend. 
		"""
		plt.figure()
		for label, data in self.lines.iteritems():
			plt.plot(self.time, numpy.array(data), label=label)
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

