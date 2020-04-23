all: plot3a.pdf plot3b.pdf

plot3a.pdf: plot3a.py data3a.txt
	python plot3a.py

plot3b.pdf: plot3b.py data3b.txt
	python plot3b.py

data3a.txt data3b: Aufgabe3.cpp
	g++ Aufgabe3.cpp -o euler
	./euler

clean:
	rm euler
