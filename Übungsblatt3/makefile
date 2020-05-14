
INCPATH = -I /usr/include/eigen3

all: Aufgabe2.txt Aufgabe1b_EV_EW.txt integration.txt

integration.txt: Aufgabe3.cpp
	g++ $(INCPATH) Aufgabe3.cpp -o Aufgabe3
	./Aufgabe3

Aufgabe2.txt: Aufgabe2.cpp
		g++ $(INCPATH) Aufgabe2.cpp -o Aufgabe2
		./Aufgabe2

Aufgabe1b_EV_EW.txt: Aufgabe1.cpp
		g++ $(INCPATH) Aufgabe1.cpp -o Aufgabe1
		./Aufgabe1


clean:
	rm -f Aufgabe1
	rm -f Aufgabe2
	rm -f Aufgabe2
