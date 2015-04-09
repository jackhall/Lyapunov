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

from system_manipulation import *
from simulation_utils import *
from common_subsystems import *
from signals import *
from solvers import *
import demo

