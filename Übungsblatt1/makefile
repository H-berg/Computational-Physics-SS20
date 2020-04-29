all: data2c.txt plot2d.pdf

INCPATH = -I /usr/include/eigen

data2c.txt: Aufgabe2.cpp
	g++ $(INCPATH) Aufgabe2.cpp -o Aufgabe2
	./Aufgabe2

plot2d.pdf: Aufgabe2d.py

Aufgabe2d.py: data2a.txt data2c.txt
	python Aufgabe2d.py

clean:
	rm Aufgabe1
