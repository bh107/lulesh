#CXX=clang++-3.5
CXX=g++
DEBUG=
#CXXFLAGS=-Wall -Wextra -O3 -std=c++11 $(DEBUG)
CXXFLAGS=-Wall -Wextra -O3 -std=c++11 -fopenmp $(DEBUG)
LIBS=-lm -L$(HOME)/.local/lib -lbh
EXTRAS=
INCLUDE+=-I$(HOME)/.local/include -I$(HOME)/.local/include/bohrium

all: lulesh_bohrium lulesh_bohrium_v2

lulesh_bohrium: src/lulesh_bohrium.cpp
	mkdir -p bin
	$(CXX) $(UTILS) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) $(EXTRAS) -o bin/$@ 

lulesh_bohrium_v2: src/lulesh_bohrium_v2.cpp
	mkdir -p bin
	$(CXX) $(UTILS) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) $(EXTRAS) -o bin/$@ 

clean:
	rm -rf *.o *~
	find bin/ -type f -name "*" -exec rm  {} \;
