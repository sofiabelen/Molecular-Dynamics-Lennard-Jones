CPP=g++
FLAGS=-O2 -Wall 
EGGFLAGS=-I/usr/local/include  -L/usr/local/lib   -leggx -lX11 -lm
COMPILER=egg

all: simulation

simulation: main.o vector.o system.o output.o gui.o
	$(CPP) $(FLAGS) main.o vector.o system.o output.o gui.o -o simulation $(EGGFLAGS)

main.o: main.cpp
	$(CPP) $(FLAGS) -c main.cpp $(EGGFLAGS)

vector.o: vector.cpp
	$(CPP) $(FLAGS) -c vector.cpp

system.o: system.cpp
	$(CPP) $(FLAGS) -c system.cpp

output.o: output.cpp
	$(CPP) $(FLAGS) -c output.cpp

gui.o: gui.cpp
	$(CPP) $(FLAGS) -c gui.cpp $(EGGFLAGS)

clean:
	rm *.o simulation
