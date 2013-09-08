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
import itertools as it
import operator
import numpy
import inspect
import matplotlib.pyplot as plt

#################
# System manipulation
State = collections.namedtuple('State', 't, x')


class ParallelSystems(object):
    """
    A container/manager for subsystems, itself a system.

    ParallelSystems acts as an aggregate of other objects which implement 
    some or all of the system concept. All subsystems are expected to have a 
    `state` tuple attribute as in `(t, x)`, but the 'x' member may be an empty 
    tuple if the subsystem has no state. Such subsystems do not need to be 
    callable, but all others are. Each subsystem object must still map states 
    and derivatives one-to-one.

    Events are represented by an instance of ParallelEvents.

    You may find it convenient to implement most supersystems by deriving from 
    ParallelEvents and simply extending `__init__`.  If so, make sure to call 
    the original `ParallelSystems.__init__` during initialization. 
    """
    def __init__(self, sys_list):
        """
        Initialize with an iterable of the subsystems you wish to manage. 
        
        Usage: ParallelSystems(subsystem_sequence)
        
        The order in which subsystems are given will determine in which order 
        state is accessed and the subsystems are called. ParallelSystems will 
        print a warning if the subsystems are not time-synchronized. In the 
        case of conflict, it will always take the current time of the first 
        subsystem to be the correct value. If any of the systems have state 
        but are not callable, `__init__` will raise a NotImplementedError.
        """
        self._subsystems = list(sys_list) #use an OrderedDict?
        self._dof = [len(sys.state[1]) for sys in sys_list]
        self._time = sys_list[0].state[0]
        self.events = ParallelEvents(sys_list)
        if not self.are_synchronized():
            print "Subsystem times aren't synchronized yet."
        #Check to make sure that any system that has state is callable.
        for dof, can_call in zip(self._dof, map(callable, sys_list)):
            if dof > 0 and not can_call:
                raise NotImplementedError("Systems with states " 
                        + "should be callable")

    def __call__(self):
        """
        Call each callable subsystem in turn and concatenate the results.
        """
        call_iter = it.ifilter(callable, self._subsystems)
        return tuple(
                 it.chain.from_iterable(
                   sys() for sys in call_iter))

    @property
    def state(self):
        """
        Collect state tuples from each subsystem in order and concatenate them. 
        For the time element it will return an internal time attribute as 
        initialized from the first subsystem. 
        """
        return State(self._time, 
                     tuple(it.chain.from_iterable(
                             sys.state[1] for sys in self._subsystems)))

    @state.setter
    def state(self, t_x):
        """
        Slice the given state tuple according to how many states each 
        subsystem had when `ParallelSystems` was initialized. The given time 
        will be written to each subsystem. Thus setting `state` will always 
        result in synchronized subsystems - time is not relative. 
        """
        self._time, x = t_x
        a = 0
        if len(self._dof) != len(x):
            raise IndexError("Length of state tuples don't match")
        for dof, sys in it.izip(self._dof, self._subsystems):
            b = a + dof 
            sys.state = self._time, x[a:b]
            a = b

    def are_synchronized(self, index=None):
        """
        Return false if any subsystems do not match the instance's internal
        record of time. Specifying an index will check only the indexed 
        subsystem.

        Usage: sync_bool = systems.are_synchronized()
               synv_bool = systems.are_synchronized(index)
        """
        if index is None:
            return not any(sys.state[0] != self._time 
                           for sys in self._subsystems)
        else:
            return self._subsystems[index].state[0] == self._time

    def synchronize(self, time=None):
        """
        Set the time element of `state` for each subsystem to 'time'. If no
        time is given, the instance will use its internal time data.

        Usage: systems.synchronize()
               systems.synchronize(time_float)
        """
        time = time or self._time
        for sys in self._subsystems:
            sys.state = time, sys.state[1]

    def add_subsystem(self, index, new_system):
        """
        Starts managing a new subsystem, inserted into the update order
        at index. All subsequent subsystems are shifted back. If the new
        subsystem's current time does not match the instance's, the subsystem
        will be synchronized and this method will print a warning. If the
        subsystem has states but is not callable, `add_subsystem` will raise
        a NotImplementedError.

        Usage: systems.add_subsystem(index, new_system)
        """
        if len(new_system.state[1]) > 0 and not callable(new_system):
            raise NotImplementedError("Systems with states should be callable")
        self._subsystems.insert(index, new_system)
        self._dof.insert(index, len(new_system.state[1]))
        if not self.are_synchronized(index):
            print "warning - forced to synchronize new subsystem"
            self.synchronize(index)

    def remove_subsystem(self, index):
        """
        Stop managing the subsystem at `index` in the update order. All
        subsequent subsystems are moved up.

        Usage: systems.remove_subsystem(index)
        """
        self._subsystems.pop(index)
        self._dof.pop(index)


class ParallelEvents(object):
    """ 
    An iterable that chains together subsystem `events` iterables. 

    `ParallelEvents` stores references to the subsystem objects, not whatever 
    object is bound to their `events` attributes. This makes updating 
    subsystem event iterables much easier, but note that this class will not 
    recheck for an `events` attribute after being initialized. If you want to 
    rerun this check, create a new instance.

    Remember that each event should be callable with no arguments.
    """
    def __init__(self, sys_list):
        """
        The new instance will check each subsystem for an `events` attribute
        and keep a reference of subsystems that have one. Order of access is
        preserved.

        Usage: events = ParallelEvents(subsystem_sequence)
        """
        def has_events(sys):
            try:
                sys.events
                return True
            except AttributeError:
                return False
        self._subsystems = filter(has_events, sys_list)

    def __len__(self):
        """ 
        Sum the lengths of constituent event iterables. Will throw if any of 
        them lack `__len__`.
        """
        return sum(len(sys.events) for sys in self._subsystems)

    def __iter__(self):
        """
        Yield at each event function in order.
        """
        for system in self._subsystems:
            for event in system.events:
                yield event


def state_property(tname='_lyapunov__t', xname='_lyapunov__x'):
    """ 
    Returns a property that acts like a `State` namedtuple. Regarless of what
    kind of tuple is assigned to it, this property will always yield a suitable
    namedtuple when accessed. The zeroth element of `State` can be referred to
    as `t` and the first as `x`.

    The elements are actually stored in attributes of the system instance. The 
    name of each attribute is mangled by default, but you can specify you own
    names for convenience. If you do choose your own names, treat those 
    attributes as read-only.

    Usage: 
        class System(object):
            state = state_property(tname='_lyapunov__t', xname='_lyapunov__x')
    """
    def fget(obj):
        return State(obj.__dict__[tname], obj.__dict__[xname])
    def fset(obj, t_x):
        obj.__dict__[tname], obj.__dict__[xname] = t_x
    def fdel(obj):
        del obj.__dict__[tname], obj.__dict__[xname]
    doc = "Acts like a StateTuple."
    return property(fget, fset, fdel, doc)


def systemfunctor(system_class):
    """
    Decorator that adds defaulted time and state arguments to a system class. 
    This should make it compatible with both Lyapunov and SciPy.
    """
    argspec = inspect.getargspec(system_class.__call__)
    if len(argspec.args) != 1:
        raise NotImplementedError("System __call__ has arity > 0")
    original_call = system_class.__call__
    def new_call(self, t=None, x=None):
        if t is not None and x is not None:
            self.state = t, x
        return original_call(self)
    system_class.__call__ = new_call
    return system_class

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
        self.events = [lambda: step_time - self.time]
        self.update()
    
    state = state_property("time", "_lyapunov__x")

    def update(self):
        if self.events[0]() > 0:
            self.value = self.initial
        else:
            self.value = self.final
    

class SquareWave(object):
    """
    Generates a square wave. The value depends on '[SquareWave].time'.
    """
    def __init__(self, period=1.0, y_lower=-1.0, y_upper=1.0):
        """
        Initialize with the period, the trough value, and the peak value.

        Usage: SquareWave(period=1.0, y_lower=-1.0, y_upper=1.0)
        """
        self.epsilon = period / 100
        self.period, self.lower, self.upper = period, y_lower, y_upper
        self.state = 0.0, ()
        self.events = [lambda: math.sin(2 * math.pi * 
                                        (self.time + self.epsilon) / period)]
        self.update()

    state = state_property("time", "_lyapunov__x")

    def update(self):
        """ Updates value - y_upper if event()>0, y_lower if not. """
        if self.events[0]() > 0:
            self.value = self.upper
        else:
            self.value = self.lower


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

#################
# Simulation Utilties
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


def check_NaN(system):
    """ 
    Check to make sure system states are still numbers. 
    Raises an exception otherwise.
    """
    if any( map(math.isnan, system.state[1]) ):
        raise ArithmeticError("System state is NaN!")

#################

