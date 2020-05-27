CXXFLAGS = -O3 -ffast-math -march=native
INCPATH = -I /usr/include/eigen3

all: A2.txt

A2.txt: Aufgabe02.cpp
	g++ $(INCPATH) $(CXXFLAGS) -std=c++11 Aufgabe02.cpp -o Aufgabe2
	./Aufgabe2

clean:
	rm -f Aufgabe1
	rm -f Aufgabe2
