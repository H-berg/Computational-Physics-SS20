CXXFLAGS = -O3 -ffast-math -march=native
INCPATH = -I /usr/include/eigen3

all: plot2b.pdf

plot2b.pdf: data2b.txt Aufgabe2.py
	python Aufgabe2.py

data2b.txt: Aufgabe2.cpp profiler.cpp profiler.h
		g++ $(CXXFLAGS) $(INCPATH) Aufgabe2.cpp profiler.cpp -o Aufgabe2
		./Aufgabe2

clean:
	rm -f ./profilerExample

# .PHONY gibt an, dass clean auch funktioniert, wenn eine Datei namens 'clean' existiert
.PHONY: clean
