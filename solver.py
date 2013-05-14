import math
import numpy
import matplotlib.pyplot as plt
import lyapunov

class System:
	"""A state machine archetype representing a dynamic system."""
	def __init__(self, state0):
		self.state = state0
		self._num_states = len(state0)

	def __len__(self):
		return self._num_states

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self, x):
		"""Extend this to compute output signals (self._output)!
			Call this version with System.state.fset(self, x)."""
		if len(x) != self._num_states:
			raise RuntimeError(str(type(self)) + ' has ' + str(self._num_states)
							   + ' states, not ' + str(len(x)) + '.')
		self._state = x

	def __call__(self, signal, disturbance=None):
		"""Returns state derivatives in a list, as required by vode.
		   Should not change the state of the system."""
		if self.__call__.__func__ is System.__call__.__func__:
			raise NotImplementedError(str(type(self)) 
							   		  + ' does not override __call__.')


class DemoNoEvents(System):
	"""Mass spring damper system."""
	def __init__(self):
		System.__init__(self, [1.0, 1.0])

	def __call__(self):
		x, v = self.state
		return [v, -(v + x)]

	def plot(self):
		plt.figure()
		try:
			plt.plot(self.x_out[:,0], self.x_out[:,1])
			plt.show()
		except:
			plt.close()


class DemoEvents(System):
	"""Double integrator with linear switching mode."""
	def __init__(self):
		System.__init__(self, [1.0, 1.0])
	
	@property
	def mode(self):
		x, v = self.state
		return v < -0.5*x

	def __call__(self):
		x, v = self.state
		u = 1 if self.mode else -1
		return [v, u]

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


class ReferenceFilter(System):
	"""Differentiates a reference signal with a linear filter. The 
		'output' attribute refers to the complete reference signal."""
	def __init__(self, gains, reference0):
		"""Place poles before constructing. 'gains' is a list of
			coefficients of the characteristic equation (normalized),
			from lowest-highest order. Exclude the highest, since it 
			should be equal to one anyway."""
		#This way, there's no need to handle complex numbers.
		state0 = [reference0].append([0.0]*(len(gains) - 1))
		System.__init__(self, state0)
		self._gains = gains #(244, 117.2, 18.75) #check signs?

	def __call__(self, signal):
		"""Since the input signal is needed when the state is set,
			there's no need to pass it to __call__."""
		xndot = (sum(-q*d for (q,d) in zip(self._state, self._gains)) 
				 + self._gains[0]*signal) #computes the only nontrivial derivative
		return self._state[1:].append(xndot)



#def ode_func(time, state, system):
#	system.state = state
#	system.time = time
#	return system()


class TimeInterval:
	def __init__(self, t1, t2):
		self.lower = t1
		self.upper = t2
	@property
	def length(self):
		return self.upper - self.lower
	@property
	def midpoint(self):
		return (self.upper + self.lower) / 2.0


def norm(x):
	result = 0.0
	for i in x:
		result += i**2
	return math.sqrt(result)


class Solver:
	def __init__(self, system, events=False, slide=True, min_ratio=.01): 
		self.system = system
		self.stepper = lyapunov.Stepper(system)
		#Check basic requirements of a system object...
		assert hasattr(system, 'state') 	#for setting state
		try: 
			system() 	#for computing derivatives
		except: 
			raise TypeError('system must be callable without arguments')
		self.events = events is True
		self.slide = slide is True
		self._autonomous = not hasattr(system, 'time') 
		self.num_points = 100 #number of points recorded/plotted
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

	def _find_root(self, step_size):
		"""A bisection rootfinder."""
		interval = TimeInterval(self.system.time - step_size, 
								self.system.time)
		min_step_size = interval.length * self.min_ratio 
		end_x, end_t = list(self.system.state), self.system.time
		self.stepper.revert() #take one step backwards
		old_mode = self.system.mode
		while interval.length > min_step_size: #and norm(self.system.state) > 0.1:
			#print (interval.lower, interval.upper)
			self.stepper.step(interval.length / 2.0)
			assert interval.length > min_step_size*0.1
			if self.system.mode == old_mode:
				interval.lower = interval.midpoint
			else:
				interval.upper = interval.midpoint
				end_x, end_t = list(self.system.state), self.system.time
				self.stepper.revert()
		if self.system.mode == old_mode:
			self.system.state, self.system.time = end_x, end_t
		

	def _slide(self, step_size, final_time):
		min_step_size = step_size * self.min_ratio
		self.stepper.step(min_step_size)
		old_mode = self.system.mode
		self.stepper.step(min_step_size)
		while self.system.mode != old_mode and self.system.time < final_time:
			old_mode = self.system.mode
			self.stepper.euler_step(min_step_size)
			self.x_out.append(self.system.state)
			self.t_out.append(self.system.time)
		return self.system.time >= final_time

	def simulate(self, final_time):
		if self._autonomous:
			self.system.time = 0.0 #means that system.time is reserved!
		step_size = (final_time - self.system.time) / self.num_points
		self.x_out = [list(self.system.state)]
		self.t_out = [self.system.time]
		if self.events is True:
			current_mode = self.system.mode
		#main solver loop
		while self.system.time < final_time:
			self.stepper.step(step_size)
			if self.events is True:
				#Detect an event - defined as a change in system.mode.
				#Currently assumes only two modes!
				if self.system.mode != current_mode:
					self._find_root(step_size)
					current_mode = self.system.mode
					if self.slide is True:
						self.stepper.step(step_size)
						if self.system.mode != current_mode:
							self.stepper.revert()
							if self._slide(step_size, final_time): break
			#Record for output.
			self.x_out.append(list(self.system.state))
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
			self.y_out.append(list(self.system.output))
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

