#!/usr/bin/python

from distutils.core import setup, Extension

#boost 1.53 or later is required
#make sure the boost install includes boost.python and boost.numeric.odeint

setup(name="lyapunov",
      version="3.0.1",
      description="a library for integrating hybrid ODEs",
      author="Jack Hall",
      author_email="jackwhall7@gmail.com",
      url="https://github.com/jackhall/Lyapunov",
      license="GNU-GPLv3",
      data_files=[("", ["GNU-GPLv3"])],
      requires=["numpy", "matplotlib"], #not sure which versions are ok
      packages=["lyapunov", "lyapunov.demo"],
      ext_modules=[Extension("lyapunov.solvers", 
                             ["lyapunov/solvers.cpp"],
                             include_dirs=["include"], 
                             libraries=["boost_python", "boost_system"],
                             extra_compile_args=["-std=c++11"],
                             )],
      )

