all: data1.txt data2c.txt plot2d.pdf

INCPATH = -I /usr/include/eigen3

data1.txt: Aufgabe1.cpp
	g++ $(INCPATH) Aufgabe1.cpp -o Aufgabe1
	./Aufgabe1
	
data2c.txt: Aufgabe2.cpp
	g++ $(INCPATH) Aufgabe2.cpp -o Aufgabe2
	./Aufgabe2

plot2d.pdf: Aufgabe2d.py

Aufgabe2d.py: data2a.txt data2c.txt
	python Aufgabe2d.py

clean:
	rm Aufgabe1
