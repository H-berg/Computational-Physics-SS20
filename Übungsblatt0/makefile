all: plota.pdf plotb.pdf

plota.pdf: plota.py eulera.txt
	python plota.py

plotb.pdf: plotb.py eulerb.txt
	python plotb.py

eulera.txt: Aufgabe3.cpp
	g++ Aufgabe3.cpp -o euler
	./euler

eulerb.txt: Aufgabe3.cpp
	g++ Aufgabe3.cpp -o euler
	./euler

clean:
	rm euler
