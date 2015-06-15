BENCH_SRC   = $(wildcard src/*.cpp)
BENCHMARK   = $(subst .cpp,, $(subst src/,,$(BENCH_SRC)))

#CXX=clang++-3.5
CXX=g++
DEBUG=
#CXXFLAGS=-Wall -Wextra -O3 -std=c++11 $(DEBUG)
CXXFLAGS=-Wall -Wextra -O3 -std=c++11 -fopenmp $(DEBUG)
LIBS=-lm -L$(HOME)/.local/lib -lbh
EXTRAS=
INCLUDE+=-I$(HOME)/.local/include -I$(HOME)/.local/include/bohrium

all: $(BENCHMARK)

$(BENCHMARK): $(BENCH_SRC)
	mkdir -p bin
	$(CXX) $(UTILS) $(CXXFLAGS) $(INCLUDE) $< $(LIBS) $(EXTRAS) -o bin/$@ 

clean:
	rm -rf *.o *~
	find bin/ -type f -name "*" -exec rm  {} \;
