INCPATH = -I /usr/include/eigen3

all: Aufgabe3.cpp Energie.pdf Energie_Equi_geordnet.pdf Energie_Equi.pdf Magnetisierung.pdf Betrag_Magnetisierung.pdf momentaufnahme_1.pdf momentaufnahme_2.pdf momentaufnahme_3.pdf

Energie.pdf Energie_Equi_geordnet.pdf Energie_Equi.pdf Magnetisierung.pdf Betrag_Magnetisierung.pdf momentaufnahme_1.pdf momentaufnahme_2.pdf momentaufnahme_3.pdf: Energie_1.txt Energie_2.txt Energie_3.txt Energie_Equi_geordnet_1.txt Energie_Equi_geordnet_2.txt Energie_Equi_geordnet_3.txt Energie_Equi_1.txt Energie_Equi_2.txt Energie_Equi_3.txt Magnetisierung_1.txt Magnetisierung_2.txt Magnetisierung_3.txt Betrag_Magnetisierung_1.txt Betrag_Magnetisierung_2.txt Betrag_Magnetisierung_3.txt spinup_1.txt spinup_2.txt spinup_3.txt spindown_1.txt spindown_2.txt spindown_3.txt simulation.py
	python simulation.py

Energie_1.txt Energie_2.txt Energie_3.txt Energie_Equi_geordnet_1.txt Energie_Equi_geordnet_2.txt Energie_Equi_geordnet_3.txt Energie_Equi_1.txt Energie_Equi_2.txt Energie_Equi_3.txt spinup_1.txt spinup_2.txt spinup_3.txt spindown_1.txt spindown_2.txt spindown_3.txt: simulation.cpp
	g++ $(INCPATH) -std=c++11 simulation.cpp -o simulation
	./simulation

Aufgabe3.cpp:
	g++ $(INCPATH) Aufgabe3.cpp -o Aufgabe3
		./Aufgabe3

clean:
	rm Aufgabe3
	rm -f ./simulation
