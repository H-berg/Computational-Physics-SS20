CXXFLAGS = -O3 -ffast-math -march=native
INCPATH = -I /usr/include/eigen3

all: Aufgabe02

a_1.cpp:
	g++ $(INCPATH) $(CXXFLAGS) -std=c++11 a_1.cpp -o a_1
	./a_1

Aufgabe02: in.txt out.txt Aufgabe02_plot.py
	python Aufgabe02_plot.py

in.txt out.txt: Aufgabe02.cpp
	g++ $(INCPATH) $(CXXFLAGS) -std=c++11 Aufgabe02.cpp -o Aufgabe02
	./Aufgabe02

clean:
	rm -f a_1
	rm -f Aufgabe02
