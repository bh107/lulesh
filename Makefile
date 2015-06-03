ROOT=/root/.local
HEADER=$(ROOT)/include/*
CPPB_INCLUDE=$(ROOT)/bridge/cpp
BIN_DIR = bin

CXX=clang++-3.5
EXTRAS+=
CXXFLAGS=-Wall -Wextra -pedantic -g -O2 -std=c++0x $(EXTRAS)

LIBS+=-L$(ROOT)/lib -lbh
INCLUDE+=-I$(ROOT)/include -I$(ROOT)/include/bohrium

all : lulesh

lulesh : src/lulesh_bohrium.cpp $(HEADER)
	$(CXX) $< -o bin/$@ -L$(ROOT)/core -I$(ROOT)/include -I$(CPPB_INCLUDE) $(INCLUDE) $(LIBS) -lbh $(LCFLAGS) $(CXXFLAGS) -lstdc++

clean :
	rm -f $(BIN_DIR)/* *.o
