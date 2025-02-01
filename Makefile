CPP=clang++
CPPFLAGS= -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp -std=c++11 -mcpu=apple-m1 -march=native -O3 

default: rt

all: rt

rt: rt.cpp Makefile
	$(CPP) $(CPPFLAGS) rt.cpp -o rt

clean:
	-rm rt
