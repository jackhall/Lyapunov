import math
import numpy
import matplotlib.pyplot as plt
import solvers


class Output(object):
	""" A descriptor for block-diagram-style outputs. """
	def __init__(self, fget):
		self.fget = fget
		self._sinks = []

	def __get__(self, obj, objtype=None):
		if obj is None:
			return self
		return self.fget(obj)

	def __delete__(self, obj):
		for sink in self._sinks:
			sink.break_link()

	def add_sink(self, sink):
		if not hasattr(sink, "break_link"):
			raise AttributeError("Sinks should implement a 'break_link' method.")
		self._sinks.append(sink)


#instance descriptor ideas:
#	assign callback to a "hidden" name in each instance
#	make descriptor a metaclass

class Input(object):
	""" A descriptor for block-diagram-stype inputs. """
	def __init(self, instance, attribute):
		self.instance = instance
		self.attribute = attribute

	def __get__(self, obj, objtype=None):
		if obj is None:
			return self
		return getattr(self.instance, self.attribute)
	
	def link_to(self, instance, attribute):
		""" Meant to be called as a method, not as a decorator. """
		self.instance = instance #sets descriptor behavior for all instances!
		self.attribute = attribute
		#self.fget = source.fget
		#source.add_sink(self)

	def break_link(self):
		self.fget = None

#write mode descriptor/decorator?


class SimpleDemo(object):
	"""Mass spring damper system."""
	def __init__(self):
		self.state = (1.0, 1.0) 

	def __len__(self):
		return 2

	def __call__(self):
		_, v = self._state
		return (v, self.a)

	u = Input()

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, x):
		self._state = x
		try:
			self.a = -x[1] - x[0] + self.u
		except NotImplementedError:
			self.a = -x[1] - x[0]

	def plot(self):
		plt.figure()
		try:
			plt.plot(self.x_out[:,0], self.x_out[:,1])
			plt.show()
		except:
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
		except:
			plt.close()


class ControlDemo(object):
	pass



class Filter(object):
	""" Differentiates a reference signal with a linear filter. The 
		'output' attribute refers to the complete reference signal. """
	def __init__(self, gains, signal0):
		""" Place poles before constructing. 'gains' is a list of
			coefficients of the characteristic equation (normalized),
			from lowest-highest order. Exclude the highest, since it 
			should be equal to one anyway. """
		#This way, there's no need to handle complex numbers.
		self._num_states = len(gains)
		self.state = (signal0,) + (0.0,)*(self._num_states - 1)
		self._xndot = 0.0
		self._gains = gains #(244, 117.2, 18.75) #check signs?

	def __len__(self):
		return self._num_states

	signal = Input()

	@Output
	def output(self):
		return self.state + (self._xndot,)

	def __call__(self):
		""" Computes the only nontrivial derivative. """
		self._xndot = (sum(-q*d for (q,d) in zip(self.state, self._gains)) 
					+ self._gains[0]*self.signal) 
		return self.state[1:] + (self._xndot,)



#def ode_func(time, state, system):
#	system.state = state
#	system.time = time
#	return system()


#class TimeInterval(object):
#	def __init__(self, t1, t2):
#		self.lower = t1
#		self.upper = t2
#	@property
#	def length(self):
#		return self.upper - self.lower
#	@property
#	def midpoint(self):
#		return (self.upper + self.lower) / 2.0


def norm(x):
	result = 0.0
	for i in x:
		result += i**2
	return math.sqrt(result)


class Solver(object):
	def __init__(self, system, events=False, slide=True, min_ratio=.01, 
			     points=100): 
		self.system = system
		self.stepper = solvers.Stepper(system)
		#Check basic requirements of a system object...
		assert hasattr(system, 'state') 	#for setting state
		try: 
			system() 	#for computing derivatives
		except: 
			raise TypeError('system must be callable without arguments')
		self.events = events is True
		self.slide = slide is True
		self._autonomous = not hasattr(system, 'time') 
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

	#def _find_root(self, step_size):
	#	"""A bisection rootfinder."""
	#	interval = TimeInterval(self.system.time - step_size, 
	#							self.system.time)
	#	min_step_size = interval.length * self.min_ratio 
	#	end_x, end_t = self.system.state, self.system.time
	#	self.stepper.revert() #take one step backwards
	#	old_mode = str(self.system.mode)
	#	while interval.length > min_step_size: 
	#		#print (interval.lower, interval.upper)
	#		self.stepper.step(interval.length / 2.0)
	#		assert interval.length > min_step_size*0.1
	#		if self.system.mode == old_mode:
	#			interval.lower = interval.midpoint
	#		else:
	#			interval.upper = interval.midpoint
	#			end_x, end_t = self.system.state, self.system.time
	#			self.stepper.revert()
	#	if self.system.mode == old_mode:
	#		self.system.state, self.system.time = end_x, end_t

	def _slide(self, step_size, final_time):
		min_step_size = step_size * self.min_ratio
		self.stepper.step(min_step_size)
		old_mode = str(self.system.mode)
		self.stepper.step(min_step_size)
		#somehow the solver gets stuck in this loop for small step sizes
		while self.system.mode != old_mode and self.system.time < final_time:
			old_mode = str(self.system.mode)
			self.stepper.euler_step(min_step_size)
			if self.system.mode == old_mode:
				self.stepper.revert()
				self.stepper.euler_step(step_size)
			self.x_out.append(self.system.state)
			self.t_out.append(self.system.time)
		return self.system.time >= final_time

	def simulate(self, final_time):
		if self._autonomous:
			self.system.time = 0.0 #means that system.time is reserved!
		step_size = (final_time - self.system.time) / self.points
		self.x_out = [self.system.state]
		self.t_out = [self.system.time]
		if self.events is True:
			current_mode = str(self.system.mode) #assumes string type!
		#main solver loop
		while self.system.time < final_time:
			self.stepper.step(step_size)
			if self.events is True:
				#Detect an event - defined as a change in system.mode.
				if self.system.mode != current_mode:
					#self._find_root(step_size)
					#C++ version is now less buggy than the python version
					self.stepper.find_root(step_size, 
										   step_size*self.min_ratio);
					assert current_mode != self.system.mode
					current_mode = str(self.system.mode)
					if self.slide is True:
						self.x_out.append(self.system.state)
						self.t_out.append(self.system.time)
						self.stepper.step(step_size)
						if self.system.mode != current_mode:
							self.stepper.revert()
							if self._slide(step_size, final_time): break
			#Record for output.
			self.x_out.append(self.system.state)
			self.t_out.append(self.system.time)
		self.system.state, self.system.time = self.x_out[0], self.t_out[0]
		if self._autonomous:
		 	del self.system.time
		try:
			self.compute_output()
			return (numpy.array(self.x_out), 
					numpy.array(self.y_out), 
					numpy.array(self.t_out))
		except: 
			return numpy.array(self.x_out), numpy.array(self.t_out)
	
	def compute_output(self):
		assert hasattr(self.system, 'output')
		self.y_out = []
		for x, t in zip(self.x_out, self.t_out):
			self.system.state, self.system.time = x, t
			self.y_out.append(self.system.output)
		self.system.state, self.system.time = self.x_out[0], self.t_out[0]
		return numpy.array(self.y_out)

	
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
		
	#def plot_system(self):
	#	#t = numpy.linspace(0, self.final_time, self.points) #time vector
	#	plt.figure() #plot object
	#	final_time_list = self.final_time
	#	try:
    #			[ x for x in self.final_time]
	#	except TypeError:
	#		final_time_list = [ self.final_time for x in range(len(self.initial_state)) ]
	#	for state0, tf in zip(self.initial_state, final_time_list):
	#		states, time = self.simulate(state0, tf)
	#		#integrate.odeint(self.f, state0, t, args=(self,))
	#		x, y = states[:,0], states[:,1]
	#		if max(x) > xmax: 
	#			xmax = max(x)
	#		if min(x) < xmin: 
	#			xmin = min(x)
	#		if max(y) > ymax: 
	#			ymax = max(y)
	#		if min(y) < ymin: 
	#			ymin = min(y)
	#		plt.plot(x, y, color='blue')
	#		plt.plot(state0[0], state0[1], 'bo', markersize=7)
	#	self.phase_portrait(xmin - self.xmargins[0], xmax + self.xmargins[1], ymin - self.ymargins[0], ymax + self.ymargins[1])
	#	plt.xlabel(self.xlabel)
	#	plt.ylabel(self.ylabel)
	#	plt.title(self.title)
	#	plt.savefig(self.filename, dpi=300)
	#	plt.show()

	

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

