from system_manipulation import state_property


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

