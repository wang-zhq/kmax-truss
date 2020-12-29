MACRO=-DDEBUG_DF
CXX=g++
CXXFLAGS=-O3 -std=c++17 -fopenmp $(MACRO)

kmtruss : kmax_truss.o kmax_dep.o
	$(CXX) $(CXXFLAGS) -o kmtruss kmax_truss.o kmax_dep.o

kmax_truss.o : kmax_truss.cpp kmax_truss.hpp
	$(CXX) $(CXXFLAGS) -c kmax_truss.cpp -o kmax_truss.o
    
kmax_dep.o : kmax_dep.cpp kmax_truss.hpp
	$(CXX) $(CXXFLAGS) -c kmax_dep.cpp -o kmax_dep.o

.PHONY : clean
clean :
	rm kmtruss  *.o
