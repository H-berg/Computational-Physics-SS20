#include <iostream>
#include <fstream>
#include <math.h>
#include <Eigen/Dense>
#include <random>

using namespace std;
using namespace Eigen;

int main(){

  std::random_device rd;
  std::mt19937 generator(rd());

  // Vorschlagswahrscheinlichkeit
  uniform_real_distribution<double> V(0,1);
  auto randomNumber = V(generator);

  // Magnetfeld initialisieren
  VectorXd H = VectorXd::LinSpaced(1e4,-5,5);

  // Initialisierung der Energiediff, Übergangswahrscheinlichkeit,
  // Spin, bei mehreren -> Random, Magnetisierung: ∑Spins
  double dE, p;
  int Spin = 1, Magnetisierung = 0;

  std::ofstream file;
  file.open("A1.txt", std::ios_base::trunc);
  file << "# Magnetfeld H   Anzahl Magnetisierung \n";

  // Metropolis Algo
  for(int i = 0; i < H.size(); i++)
  {
    // Übergangswahrscheinlichkeit p bestimmen
    p = exp(-2*abs(H(i)));
    for(int j = 0; j < 1e5; j++)
    {
      // Energiedifferenz dE bestimmen
      dE = (-Spin - Spin)*H(i);

      // wenn dE < 0 dann akzeptieren wir den Spin-Flip, da energetisch günstiger
      if(dE < 0)
      {
        Magnetisierung += Spin;
        Spin = -Spin;
      }
      // falls dE > 0, vergleichen wir die Übergangswahrscheinlichkeit p
      // (dabei ist kT = 1) mit randomNumber ∈ [0,1]
      else
      {
        randomNumber = V(generator);
        // wir akzeptieren, wenn:
        if(randomNumber <= p)
        {
          Magnetisierung += Spin;
          Spin = -Spin;
        }
        // sonsten lassen wir den Spin unflipped
        else
        {
          Magnetisierung -= Spin;
        }
      }
    }
    // Magnetisierung speichern
    file << H(i) << "  " << Magnetisierung << endl;
    Magnetisierung = 0;
  }
  file.close();
return 0;
}
