CC=g++ 
CFLAGS= -Wall -std=c++11

executable: main.o MSERStates.o CVGeometryUtils.o
	$(CC) $(CFLAGS) -o executable  main.o MSERStates.o CVGeometryUtils.o `pkg-config --libs opencv` -lfftw3 -L./bin/GraphCutsOptimization/lib -lgraphcutsoptimization
main.o: 
	$(CC) $(CFLAGS) `pkg-config --cflags opencv` -I ./bin/GraphCutsOptimization/ -g -c main.cpp 
MSERStates.o:
	$(CC) $(CFLAGS) `pkg-config --cflags opencv` -g -c MSERStates.cpp 

CVGeometryUtils.o:
	$(CC) $(CFLAGS) `pkg-config --cflags opencv` -g -c CVGeometryUtils.cpp
