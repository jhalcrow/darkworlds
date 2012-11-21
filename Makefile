
CXX_CFLAGS=-O2 -std=c++0x

%.o: %.cpp
	c++ $(CXX_CFLAGS) -c $? -o $@

halomodel: halomodel2.o
	c++ $(CXX_CFLAGS) -o halomodel halomodel2.o

clean:
	rm halomodel halomodel2.o

all: halomodel
	#