CPP=g++
FLAGS=-O2 -Wall main.cpp -I/usr/local/include  -L/usr/local/lib   -leggx -lX11 -lm
COMPILER=egg

all:
	egg main.cpp -o simulation

clean:
	rm *.o simulation
