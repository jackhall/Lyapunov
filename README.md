Lyapunov
========

A library for numerically integrating nonlinear dynamical systems.
The library is released under the GNU General Public License version 3.
Documentation is in the code and online at https://github.com/jackhall/Lyapunov/wiki.

version 3.1.2, 7/15/2018

**Requirements**

* C++11 capable compiler
* Boost 1.53 or later (boost.python and boost.numeric.odeint)
* Python 2.7.5 or later (including Python 3.x)
* Numpy
* Matplotlib


**Compilation and Installation**

In the base directory, run `python setup.py install`.
On Windows, use `setup.py install` from the command prompt.
For other options, see the [standard Python docs](http://docs.python.org/2/install/index.html#install-index) on module installation.

Many systems have parallel installations of python. 
This can make it tricky to link the compiled solvers to the proper version of the boost libraries. 
If you have any problems getting that to work, please email me or create a github issue and I'll help you out. 
Knowing more about the environments in which people try to install Lyapunov will also help me refine the build process. 

**Uninstallation**

Remove the file "lyapunov" from wherever your python libraries are stored.
