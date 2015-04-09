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
import matplotlib.pyplot as plt


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

