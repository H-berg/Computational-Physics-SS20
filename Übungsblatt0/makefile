all: plot3a.pdf plot3b.pdf a_2_a.pdf a_2_b.pdf a_2_c.pdf a_2_rel_err.pdf

plot3a.pdf: plot3a.py data3a.txt
	python plot3a.py

plot3b.pdf: plot3b.py data3b.txt
	python plot3b.py

data3a.txt data3b: Aufgabe3.cpp
	g++ Aufgabe3.cpp -o euler
	./euler

a_2_a.txt a_2_b.txt a_2_c.txt: Aufgabe2.cpp
	g++ Aufgabe2.cpp -o Aufgabe2
	./Aufgabe2

a_2_a.pdf a_2_b.pdf a_2_c.pdf a_2_rel_err.pdf: a_2_a.txt a_2_b.txt a_2_c.txt plot2.py
	python plot2.py

clean:
	rm euler
