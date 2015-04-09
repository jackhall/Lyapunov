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
        original_call = system.__call__
        def new_call(self, t=None, x=None):
            if t is not None and x is not None:
                self.state = t, x
            return original_call(self)
        system.__call__ = new_call
        return system

