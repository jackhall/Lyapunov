CC = g++
PYTHON_INCLUDE_PATH = /usr/include/python2.7
PYTHON_LIB_PATH = /usr/lib64/python2.7/site-packages
CFLAGS = -std=c++11 -fPIC -I$(PYTHON_INCLUDE_PATH) -g
LFLAGS = -shared -lpython2.7 -lboost_python -lboost_system -Wl,-no-undefined -g

solvers.so : solvers.cpp 
	$(CC) $(CFLAGS) -o solvers.o -c solvers.cpp
	$(CC) $(LFLAGS) -o solvers.so solvers.o
	rm *.o

install :
	sudo cp solvers.so $(PYTHON_LIB_PATH)/.
	sudo cp lyapunov.py $(PYTHON_LIB_PATH)/.

remove :
	sudo rm $(PYTHON_LIB_PATH)/solvers.so $(PYTHON_LIB_PATH)/lyapunov.py

