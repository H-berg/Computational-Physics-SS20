CXXFLAGS = -O3 -ffast-math -march=native
INCPATH = -I /usr/include/eigen3

all: A2.txt a_1_1.jpg a_1_1_1.jpg a_1_1_2.jpg

a_1_1.jpg a_1_1_1.jpg a_1_1_2.jpg: dft_3.txt direkt_3.txt dft_4.txt direkt_4.txt fft_1.txt fft_2.txt
	python a_1.py

dft_3.txt direkt_3.txt dft_4.txt direkt_4.txt fft_1.txt fft_2.txt: a_1.cpp
	g++ $(CXXFLAGS) $(INCPATH) -std=c++11 a_1.cpp -o a_1
	./a_1

A2.txt: Aufgabe02.cpp
	g++ $(INCPATH) $(CXXFLAGS) -std=c++11 Aufgabe02.cpp -o Aufgabe2
	./Aufgabe2

clean:
	rm -f a_1
	rm -f Aufgabe2
