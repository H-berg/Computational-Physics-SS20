CXXFLAGS = -O3 -ffast-math -march=native
INCPATH = -I /usr/include/eigen3

all: A1 A2

A1: teil_a1.txt teil_a2.txt teil_c.txt Aufgabe01_plot.py | Abbildungen
	python Aufgabe01_plot.py

teil_a1.txt teil_a2.txt teil_c.txt: Aufgabe01.cpp
	g++ $(INCPATH) $(CXXFLAGS) -std=c++11 Aufgabe01.cpp -o Aufgabe1
	./Aufgabe1


A2: A2_1a.txt A2_2a.txt A2_3a.txt Aufgabe02_plot.py | Abbildungen
	python Aufgabe02_plot.py

A2_1a.txt A2_2a.txt A2_3a.txt: Aufgabe02.cpp
	g++ $(INCPATH) $(CXXFLAGS) -std=c++11 Aufgabe02.cpp -o Aufgabe2
	./Aufgabe2

Abbildungen:
		mkdir -p Abbildungen/

clean:
	rm -f Aufgabe1
	rm -f Aufgabe2
