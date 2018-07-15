#!/usr/bin/env python
"""
    Lyapunov is a toolbox for integrating dynamical systems. Integration of
    ordinary equations is done with solvers from boost.numeric.odeint, wrapped
    using boost.python.

    Instead of treating systems as functions, Lyapunov represents systems as
    objects. This not only significantly cleans up the solver interface, but it
    also gives the user more control over the simulation and encourages a cleaner
    coding style. There is a Tutorial, and the "demo/" directory contains several
    working demonstrations of the library interface.

    This reference starts with a function interview and then moves on to a
    detailed description of the library as a set of interface concepts.

    **Exposed Classes and Functions**

    Each of these stepper classes uses a different numerical solver. The low-level
    numerical routines are taken from boost.numeric.odeint (a widely available C++
    library) and wrapped using boost.python. Documentation for steppers is general
    - all steppers are used similarly, but here is an overview. The first four
    stepper classes implement fixed_step algorithms. Note: the order of a solver
    is a measure of how fast numerical error accumulates in the solution.

        euler - First-order solver useful mostly for demonstration purposes
        modified_midpoint - Second-order Runge-Kutta solver - simple but inaccurate
        runge_kutta4 - Commonly used fourth-order solver
        adams_bashforth* - Multistep algorithms - good for expensive system calls

    These next three include step size control. Each solver provides an
    approximate measure of the error accrued in each step, allowing Lyapunov to
    adjust the step size in response. They can also emulate fixed-step solutions
    without sacrificing error control, in exchange for a little speed.

        cash_karp - A fifth-order RK solver with fourth-order error estimation
        dormand_prince - Another fifth-order RK method with fourth-order error
                 estimation - has internal state
        fehlberg87 - An eighth-order RK solver with seventh-order error estimation

    These utilities should help the user adhere to the system concept without
    sacrificing convenience or coding style. In particular it should be easy to
    recursively nest subsystems.

        ParallelSystems - A container and manager for subsystems, itself a full
                          system
        ParallelEvents - A generator to chain together subsystem event iterables
        State - A namedtuple to help with the interface to system objects: 't'
                for element 0 and 'x' for element 1
        state_property - Returns a convenient property that acts like a State
        systemfunctor - A class decorator that modifies a system for compatibility
                        with SciPy

    Lyapunov provides some commonly used subsystem classes. The signal classes
    have no state and are completely time-dependent. Feel free to submit pull
    requests if you write a subsystem you like!

        Filter - A linear filter of arbitrary order
        PID - A basic Proportional Integral Differential controller
        StepSignal - Generates a step signal
        SquareWave - Generates a square wave
        SineWave - Generates a sinusoid
        ChirpSignal - Generates a sinusoid with an arbitrary instantaneous
                      frequency

    These utilities provide shortcuts for general use in simulations. Again, feel
    free to submit your own code!

        Recorder - A way to track and later plot system state and arbitrary data
        check_NaN - Function that takes a system and raises an ArithmeticError if
                    any states have been corrupted

    **System Objects**

    System objects must provide a state attribute that is a tuple (t, x) where t
    is the current system time and x is a tuple of the current system state values.
    Lyapunov uses tuples because they are immutable; if state was mutable then it
    could be changed without explicitly assigning it to the system, which would
    undermine the use of a state descriptor and potentially violate numerical
    assumptions of continuity. To avoid the need to access time and state by index
    - which would be ugly - Lyapunov provides a state_property descriptor that acts
    like a namedtuple.

    A system must also behave like an arity-zero function. The value returned from
    calling the system should be a tuple of floats corresponding to state
    derivatives; states and derivatives should be mapped one-to-one. The
    derivatives should be a continuous function of only state and time. This rule
    may be violated when an event occurs, however.

    When an event does occur, altering the system object directly is fair game.
    This is a bad idea during normal integration because such changes are
    generally discontinuous, but an event tells the stepper to expect
    discontinuity. For instance, you may alter the system's events attribute as
    you please, or change the state equations. If you wish to set a new state,
    use the step_across method; otherwise this particular type of change will be
    forgotten when solving continues. For now, do not change the number of states
    a system has. Later I may add support for resizing the state tuple.

    Specify events through a system object's events attribute - an iterable of
    arity-zero functions or functors. A discrete event occurs when one or more of
    these event functions changes sign. As with the state derivatives, the value
    of an event function should be unique given system state and time. Lyapunov's
    bisection method does not currently require event functions to be continuous,
    but this may change with the addition of a more advanced rootfinder.

    Many of these rules may seem overly restrictive, but they actually apply to
    any numerical integration algorithms you've used before. I may simply be more
    meticulous about including them in documentation. Lyapunov is designed with
    the Zen of Python in mind: "There should only be one obvious way to do it."
    State integration and events are orthogonal features, respectively
    representing continuous and discrete behavior. Using them as such will give
    you the best results.

    Many of these rules loosen their grip where subsystems are concerned. For
    instance, grouping subsystems together in an object-oriented way makes little
    sense unless they are coupled, in which case their behavior is no longer
    purely a function of their state and time. This is perfectly fine so long as
    the supersystem - the system object directly managed by the stepper - does
    follow the rules. In particular, callbacks provide a useful way of passing
    information back and forth as long as no information outside the supersystem
    is used to determine behavior. Make sure your subsystems are properly
    encapsulated. Lyapunov provides several classes that make nesting system
    objects much easier.

    **Planned Features**

        a noise signal for continuous randomness
        a state machine interface for discrete state (on top of the events interface)
        initial-value PDE solving with `numpy.ndarray` as the state array
        symbolic system representations with sympy
        implicit solvers
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
        self._subsystems = list(sys_list)  #use an OrderedDict?
        self._dof = [len(sys.state[1]) for sys in sys_list]
        self._time = sys_list[0].state[0]
        self.events = ParallelEvents(sys_list)
        if not self.are_synchronized():
            print("Subsystem times aren't synchronized yet.")

        #Check to make sure that any system that has state is callable.
        for dof, can_call in zip(self._dof, map(callable, sys_list)):
            if dof > 0 and not can_call:
                raise NotImplementedError("Systems with states should be"
                                          "callable")

    def __call__(self):
        """
        Call each callable subsystem in turn and concatenate the results.
        """
        call_iter = it.ifilter(callable, self._subsystems)
        return tuple(
            it.chain.from_iterable(sys() for sys in call_iter)
        )

    @property
    def state(self):
        """
        Collect state tuples from each subsystem in order and concatenate them.
        For the time element it will return an internal time attribute as
        initialized from the first subsystem.
        """
        x = it.chain.from_iterable(sys.state[1] for sys in self._subsystems)
        return State(self._time, tuple(x))

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
        if sum(self._dof) != len(x):
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
            print("warning - forced to synchronize new subsystem")
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


def systemfunctor(system):
    """
    Decorator that adds defaulted time and state arguments to a system class
    or wraps a system function in a system class.
    This should make it compatible with both Lyapunov and SciPy.

    Like any decorator, systemfunctor may also be used as a function that
    takes and returns a class for function. If the given class does not have
    `__call__` defined to take no arguments, it will raise a
    NotImplementedError.
    """
    if inspect.isfunction(system):
        class system_class(object):
            def __init__(self, *parameters):
                self.parameters = parameters

            state = state_property(xname="x", tname="t")

            def __call__(self, t=None, x=None):
                if t is not None and x is not None:
                    self.state = t, x
                return system(t, x, *self.parameters)

        return system_class

    else:
        argspec = inspect.getargspec(system.__call__)
        if len(argspec.args) != 1:
            raise NotImplementedError("System __call__ has arity > 0")

        def new_call(self, t=None, x=None):
            if t is not None and x is not None:
                self.state = t, x
            return system.__call__(self)

        system.__call__ = new_call
        return system

#################

#################
# Controllers
class PID(object):
    """
    A SISO PID controller. I wrote this class mostly to demonstrate how
    controller states are properly handled. The reference `r` is drawn from
    another subsystem (maybe one of the signal objects), and the output `y`
    should be drawn from the plant object. Neither of these callbacks has a
    default value; they must be set after initialization.
    """
    def __init__(self, Kp=1.0, Ki=0.0, Kd=0.0):
        """
        Create with linear gains for the proportional, integral, and
        derivative terms, in that order. The two callback functions `r` and
        `y` are left as None, so remember to set these using another object
        before simulating. This may seem strange, but it's natural style for
        Lyapunov because otherwise you run into chicken-and-egg problems with
        your subsystems.

        Usage: controller = PID([Kp=1.0][, Ki=0.0][, Kd=0.0])
        """
        self.Kp, self.Ki, self.Kd = Kp, Ki, Kd
        self.r = self.y = None
        self.state = 0.0, (0.0,) #integral term

    state = state_property(xname="_state")

    def __call__(self):
        """
        Since the sole state is the value of the integral term, the state
        derivative is simply the error.
        """
        x, _ = self.y()
        error = self.r() - x
        return (error,)

    def u(self):
        """
        Returns the control effort for the current error, error derivative,
        and error integral. The integral is a controller state, and the other
        terms are drawn from the `y` callback. Do not call until the callbacks
        are set.
        """
        #catch Nonetype exceptions for y and r!
        xi = self.y() #y returns output derivative
        x, v = xi[0], xi[1] #output and first derivative
        error = self.r() - x
        return self.Kp*error + self.Ki*self.state.x[0] - self.Kd*v

#################

#################
# Input Signals
class StepSignal(object):
    """
    Generates a step signal. The step itself is an event that the
    stepper will detect. When the event occurs, call `update` after
    calling `stepper.step_across`. It's important that the StepSignal
    instance be on the right side of the boundary. See `update` for
    more information.

    `value` is a data attribute set to the current value of the signal.
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
        """
        Sets `value` according to the current time. Call after
        stepping through the step event boundary.
        """
        if self.events[0]() > 0:
            self.value = self.initial
        else:
            self.value = self.final


class SquareWave(object):
    """
    Generates a square wave. The value depends on '[SquareWave].time'.
    Generates a square wave. The steps themselves are events that the
    stepper will detect. When an event occurs, call `update` after
    calling `stepper.step_across`. It's important that the SquareWave
    instance be on the right side of the boundary. See `update` for
    more information.

    `value` is a data attribute set to the current value of the signal.
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
        """
        Sets `value` according to the current time. Call after
        stepping through the step event boundary.
        """
        if self.events[0]() > 0:
            self.value = self.upper
        else:
            self.value = self.lower


class SineWave(object):
    """
    Generates a sinusoid. Its value depends on the current time.
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
    Generates a sinusoid with a time-varying frequency.
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
    A linear filter of arbitrary order. It also provides derivative estimates
    for otherwise nondifferentiable inputs. The `output` function returns to
    the complete signal, including derivatives. The signal is drawn from
    another subsystem with the `signal` callback. Make sure to set this before
    using the filter.
    """
    def __init__(self, gains):
        """
        Create with an iterable over the coefficients of the filter's
        desired characteristic equation. The coefficients should be
        normalized; exclude the highest-order coefficient. It is taken to be
        one. Order the coefficients from lowest-to-highest order. If you
        know the poles you want, you can simply convolve them into the
        characteristic equation.
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
        self._xndot = sum(-q*d for (q,d) in zip(self._state, self._gains))

    def output(self):
        """
        Returns the current filtered signal and all its derivatives. Meant for
        use as a callback.
        """
        return self._state + (self._xndot + self._gains[0]*self.signal(),)

    def __call__(self):
        return self._state[1:] + (self._xndot + self._gains[0]*self.signal(),)

#################
# Simulation Utilties
class Recorder(object):
    """
    Records system data during a simulation through callbacks.

    Use Recorder to record arbitrary system data during a simulation and
    plot it afterwards. Recorder will also record a full time and state
    history in it `x` and `t` data members.

    During simulation, call `log` to record a data point, optionally including
    any events which have occurred.

    After simulation, call either 'time_response' or 'phase_portrait' to
    plot. The former will plot all recorded quantities against time. The
    latter should be called with two labels, which tell it which two
    variables to plot against each other.

    Calling 'clear' will preserve labels and callbacks, but delete all data
    specific to a given simulation run.
    """
    def __init__(self, system, labels={}):
        """
        Create a Recorder from a dict mapping variable labels to
        callback functions. The callback functions should be callable
        with no argments and return a scalar numeric.

        Usage: record = Recorder(system)
               record = Recorder(system, {label_string: labeled_function,})
        """
        self.system = system
        self.labels = labels
        self.events = {}
        self.lines = {label: [] for label in labels.keys()}
        self.x = []
        self.t = []

    def log(self, events=[]):
        """
        Call to record system variables and callback values at the
        current state and time. If events are specified, then they are
        stored associatively with the current state and time.

        Usage: record.log()
               record.log(events_list)
        """
        if events:
            self.events[len(self.t)] = events
        self.t.append(self.system.state[0])
        self.x.append(self.system.state[1])
        for label, f in self.labels.iteritems():
            self.lines[label].append(f())

    def clear(self):
        """
        Call with no arguments to erase all saved records of a simulation.
        The callback functions and the labels themselves remain, but time,
        state, events, and labeled data is deleted.

        Usage: record.clear()
        """
        self.x = []
        self.t = []
        self.events = {}
        for label in self.lines:
            self.lines[label] = []

    def time_response(self, labels=set()):
        """
        Plots all saved callack results against time. String labels are given
        in the legend. Each line will be a different color. Events are shown as
        vertical lines at the time each event or list of events occurred. Later
        I may add some facility for labeling events. I also plan to allow the
        plotting of arbitrary state variables.

        Usage: record.time_response()
        """
        plt.figure()
        for label, data in self.lines.iteritems():
            plt.plot(self.t, numpy.array(data), label=label)
        for i in self.events.keys():
            plt.axvline(self.t[i], color='r')
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

        Events are plotted as little red circles around each event point. As
        in `time_response`, there is no way yet to label these.

        Later I may add the ability to plot arbitrary state variables against
        each other. In this case I may add the underlying vector field.

        Usage: record.phase_portrait(x_label_string, y_label_string)
        """
        plt.figure()
        plt.plot(numpy.array(self.lines[xlabel]),
                 numpy.array(self.lines[ylabel]))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        for i in self.events.keys():
            plt.plot(self.lines[xlabel][i], self.lines[ylabel][i], 'ro')
        plt.show()


def check_NaN(system):
    """
    Check to make sure system states are still numbers.
    Raises an ArithmeticError otherwise.

    Usage: check_NaN(system)
    """
    if any( map(math.isnan, system.state[1]) ):
        raise ArithmeticError("System state is NaN!")

#################
