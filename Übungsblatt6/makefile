INCPATH = -I /usr/include/eigen3

all: err.pdf plotgrad_a.pdf plotgrad_b.pdf

err.pdf plotgrad_a.pdf plotgrad_b.pdf: Plot.py data_1.txt data_2.txt
	python Plot.py

data_1.txt data_2.txt: Aufgabe1.cpp
	g++ $(INCPATH) Aufgabe1.cpp -o Aufgabe1
		./Aufgabe1

clean:
	rm Aufgabe1
