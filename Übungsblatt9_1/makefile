CXXFLAGS = -O3 -ffast-math -march=native
INCPATH = -I /usr/include/eigen3

all: A1

A1: A1.txt Aufgabe01_plot.py | Abbildungen
	python Aufgabe01_plot.py

A1.txt: Aufgabe01.cpp
	g++ $(INCPATH) $(CXXFLAGS) -std=c++11 Aufgabe01.cpp -o Aufgabe1
	./Aufgabe1


Abbildungen:
		mkdir -p Abbildungen/

clean:
	rm -f Aufgabe1
