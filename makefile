all: plota.pdf plotb.pdf

plota.pdf: plota.py eulera.txt
	python plota.py

plotb.pdf: plotb.py eulerb.txt
	python plotb.py

eulera.txt: main.cpp
	g++ main.cpp -o euler
	./euler

eulerb.txt: main.cpp
	g++ main.cpp -o euler
	./euler

clean:
	rm euler
