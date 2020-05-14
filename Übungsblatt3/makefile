
INCPATH = -I /usr/include/eigen3

all: Aufgabe2.txt Aufgabe1b_EV_EW.txt

Aufgabe2.txt: Aufgabe2.cpp
		g++ $(INCPATH) Aufgabe2.cpp -o Aufgabe2
		./Aufgabe2

Aufgabe1b_EV_EW.txt: Aufgabe1.cpp
		g++ $(INCPATH) Aufgabe1.cpp -o Aufgabe1
		./Aufgabe1


clean:
	rm Aufgabe2
