CXXFLAGS = -O3 -ffast-math -march=native
INCPATH = -I /usr/include/eigen3

all: aufgabe01 plot2b.pdf

aufgabe01: bild_rang_0.txt bild_rang_10.txt	bild_rang_20.txt bild_rang_50.txt
		python aufgabe01.py

bild_rang_0.txt bild_rang_10.txt	bild_rang_20.txt bild_rang_50.txt: aufgabe01.cpp
		g++ aufgabe01.cpp -o aufgabe01 $(INCPATH)
		./aufgabe01

plot2b.pdf: data2b.txt Aufgabe2.py
	python Aufgabe2.py

data2b.txt: Aufgabe2.cpp profiler.cpp profiler.h
		g++ $(CXXFLAGS) $(INCPATH) Aufgabe2.cpp profiler.cpp -o Aufgabe2
		./Aufgabe2

clean:
	rm -f ./profilerExample

# .PHONY gibt an, dass clean auch funktioniert, wenn eine Datei namens 'clean' existiert
.PHONY: clean
