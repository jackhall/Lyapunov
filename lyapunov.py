import math
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
		except TypeError: 
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


class PID(object):
	def __init__(self, Kp=1.0, Ki=0.0, Kd=0.0):
		self.Kp, self.Ki, self.Kd = Kp, Ki, Kd
		self.r = self.y = None
		self.state = (0.0,) #integral term
	
	def __call__(self):
		x, v = self.y()
		error = self.r() - x
		self._force = self.Kp*error + self.Ki*self.state - self.Kd*v
		return (error,)

	def u(self):
		return self._force


class SubsystemDemo(object):
	def __init__(self):
		self.plant = SimpleDemo()
		self.control = PID()
		self.control.y = lambda : self.plant.state
		self.plant.u = self.control.u
		self.control.r = self.reference
		self.time = 0.0

	@property
	def state(self):
		return self.plant.state + self.control.state

	@state.setter
	def state(self, x):
		self.plant.state = x[:-1]
		self.control.state = x[-1]

	def reference(self):
		return 0.0 if self.time < 2.0 else 1.0

	def __len__(self):
		return 3

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
		self._xndot = ( sum(-q*d for (q,d) in zip(self._state, self._gains)) 
					  + self._gains[0]*self.signal()) 

	def __len__(self):
		return self._num_states

	def output(self):
		return self._state + (self._xndot,)

	def __call__(self):
		""" Computes the only nontrivial derivative. """
		return self._state[1:] + (self._xndot,)



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
					self.stepper.find_root(step_size, 
										   step_size*self.min_ratio);
					current_mode = str(self.system.mode)
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
		except AttributeError: 
			return numpy.array(self.x_out), numpy.array(self.t_out)
	
	def compute_output(self):
		self.y_out = []
		for x, t in zip(self.x_out, self.t_out):
			self.system.state, self.system.time = x, t
			self.y_out.append(self.system.output)
		self.system.state, self.system.time = self.x_out[0], self.t_out[0]

	
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

