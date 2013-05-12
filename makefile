CC = g++
CFLAGS = -std=c++11
LFLAGS = -shared -Wl,-no-undefined
PFLAGS = -lpython2.7 -lboost_python -lboost_system

stepper.so : stepper.cpp
	$(CC) $(CFLAGS) $(PATHS) -fPIC -I/usr/include/python2.7 -o stepper.o -c stepper.cpp
	$(CC) $(LFLAGS) $(PFLAGS) -o stepper.so stepper.o
	rm *.o

#build/libbenoit.so : build/IndexBase.o
#	$(CC) $(LFLAGS) -lboost_system -Wl,-soname,libbenoit.so.1 -o build/libbenoit.so.1.0 build/IndexBase.o; \
	cd build; \
	ln -s libbenoit.so.1.0 libbenoit.so.1; \
	ln -s libbenoit.so.1 libbenoit.so

#build/IndexBase.o : src/IndexBase.cpp src/IndexBase.h src/Singleton.h
#	$(CC) $(CFLAGS) $(PATHS) -fPIC -c src/IndexBase.cpp; \
	mv IndexBase.o build/.

#install : 
#	sudo mv build/libbenoit.so /usr/lib64/.; \
	sudo mv build/libbenoit.so.1.0 /usr/lib64/.; \
	sudo ldconfig
	
#clean :
#	rm build/*.o
	
#remove :
#	rm build/*; \
	sudo rm /usr/lib64/libbenoit.*; \
	sudo ldconfig
	
