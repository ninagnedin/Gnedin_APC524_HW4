.SUFFIXES: .c .o
CC=g++
CXXFLAGS = -g -Wall
LIB=-lm


OBJ_SERIAL = heat_serial.o
OBJ_OMP = heat_omp.o
OBJ_MPI = heat_mpi.o

all: heat_serial heat_omp heat_mpi

#---------------------  targets  ---------------------------------------

heat_serial: $(OBJ_SERIAL) 
	$(CXX) -o $@ $^

heat_omp: $(OBJ_OMP) 
	$(CXX) -o $@ $^

heat_mpi: $(OBJ_MPI) 
	$(CXX) -o $@ $^

#-----------------------------------------------------------------------
clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend