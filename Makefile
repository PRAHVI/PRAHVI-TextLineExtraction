CC=g++ 
CFLAGS= -Wall -std=c++11

executable: main.o mserstates.o 
	$(CC) $(CFLAGS) -o executable `pkg-config --libs opencv` -lfftw3 -L./bin/GraphCutsOptimization/lib -lgraphcutsoptimization main.o mserstates.o 
main.o:
	$(CC) $(CFLAGS) `pkg-config --cflags opencv` -I./bin/GraphCutsOptimization/include -c main.cpp

mserstates.o:
	$(CC) $(CFLAGS) `pkg-config --cflags opencv` -c MSERStates.cpp


