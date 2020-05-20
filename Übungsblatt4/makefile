INCPATH = -I /usr/include/eigen3

all: Aufgabe02

Aufgabe02: in.txt out.txt Aufgabe02_plot.py
	python Aufgabe02_plot.py

in.txt out.txt: Aufgabe02.cpp
	g++ $(INCPATH) Aufgabe02.cpp -o Aufgabe02
	./Aufgabe02

clean:
	rm -f Aufgabe1
	rm -f Aufgabe2
	rm -f Aufgabe2
