#!/usr/bin/python

from distutils.core import setup, Extension

# boost 1.53 or later is required
# make sure the boost install includes boost.python and boost.numeric.odeint

solvers = Extension("lyapunov.solvers",
                    sources=["lyapunov/solvers.cpp"],
                    include_dirs=["include"],
                    libraries=["boost_python", "boost_system"],
                    extra_compile_args=["-std=c++11"],
                    )


setup(name="lyapunov",
      version="3.1.2",
      description="a library for integrating hybrid ODEs",
      author="Jack Hall",
      author_email="jackwhall7@gmail.com",
      url="https://github.com/jackhall/Lyapunov",
      license="BSD-3-Clause",
      data_files=[("", ["BSD-3-Clause.lic", "README.md"])],
      requires=["numpy", "matplotlib"],  # not sure which versions are ok
      packages=["lyapunov", "lyapunov.demo"],
      ext_modules=[solvers],
      )
