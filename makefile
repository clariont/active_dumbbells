all: main

CFLAGS= -O3

clean:
	-rm *.o brown
        
main: genarray.h
	g++ $(CFLAGS) -o brown brownianmd.cpp -lm -lgsl -lgslcblas
