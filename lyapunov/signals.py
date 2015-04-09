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
from system_manipulation import state_property


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

