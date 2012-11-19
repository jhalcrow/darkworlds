
CXX_CFLAGS=-O2 -std=c++0x
CFLAGS=-O2 -std=c++0x

%.o: %.cpp
	c++ $(CXX_CFLAGS) -c $? -o $@

halomodel: halomodel2.o
	c++ -O2 -std=c++0x -o halomodel halomodel2.o

clean:
	rm halomodel halomodel2.o

all: halomodel
	#