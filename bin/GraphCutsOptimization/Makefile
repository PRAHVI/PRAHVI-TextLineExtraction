CC=g++
CFLAGS= -Wall 

libgraphcutsoptimization.a: GCoptimization.o LinkedBlockList.o
	ar -cvq libgraphcutsoptimization.a GCoptimization.o LinkedBlockList.o
	mkdir -p lib
	mv libgraphcutsoptimization.a lib
	mkdir -p include
	cp GCoptimization.h include

GCoptimization.o: 
	$(CC) $(CFLAGS) -c GCoptimization.cpp

LinkedBlockList.o:
	$(CC) $(CFLAGS) -c LinkedBlockList.cpp
