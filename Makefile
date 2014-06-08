CXX = g++
CXXFLAGS ?= -O3 -ansi -Wall -pedantic #-pg  -O3 #-fopenmp -Winline
OBJECTS = main.o 
LIBS = 

.PHONY: all clean

all: lbm

# generic compilation rule
%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $<

#how to link
lbm: ${OBJECTS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJECTS} ${LIBS}

clean:
	rm -f *.o *~ *.vtk lbm
