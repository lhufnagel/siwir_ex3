CXX = g++
CXXFLAGS ?= -O3 -ansi -Wall -pedantic -fopenmp 
OBJECTS = io.o main.o 
LIBS = -lnetpbm 

.PHONY: all clean

all: lbm

# generic compilation rule
%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $<

lbm: ${OBJECTS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJECTS} ${LIBS}

clean:
	rm -f *.o *~ *.vtk lbm
