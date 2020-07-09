#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <random>
#include <math.h>

using namespace std;
using namespace Eigen;

//Wichtig: Mit dieser Größe des Modells und der Anzahl an Sweeps braucht das Programm sehr lange.
//Soll es zum Testen schneller durchlaufen, kann die Größe am Anfang von main() bei N und L verkleinert werden.
//L Muss dabei die Wurzel aus N sein und es ist wichtig, dass S nicht kleiner als N ist.
//Entsprechend muss dann noch im Python Programm simulation.py die Anzahl an Sweeps angepasst werden
//(sofern dieser verändert wurde) damit die Plots erstellt werden können.

std::random_device rd;
std::mt19937 generator(rd());

void initialisierung(MatrixXd &r, const int N, const int L, int ordnung){ //dritte Zeile gibt die Spinausrichtung an
  int count = 0;
  for (int n = 0; n < sqrt(N); n++){
    for (int m = 0; m < sqrt(N); m++){
      r(0, count) = n;
      r(1, count) = m;
      if(ordnung == 1){ //Anfangsbedingung: völlig geordnete Spins
        r(2, count) = 1;
      } else { //Anfangsbedingung: zufällige Spins
        uniform_real_distribution<double> distribution_3(0, 1000);
        double p = distribution_3(generator);
        if(p < 500){
          r(2, count) = -1;
        } else{
          r(2, count) = 1;
        }
      }
      count++;
    }
  }
}

double metropolis(double spin, VectorXd naechsteNachbarn, VectorXd &energie, double position, double kT){
  double energiedifferenz = 0, e1 = 0, e2 = 0;

  for(int i = 0; i < 4; i++){ //mit Spinflip
    e2 -= (-spin)*naechsteNachbarn(i);
  }
  for(int i = 0; i < 4; i++){ //ohne Spinflip
    e1 -= (spin)*naechsteNachbarn(i);
  }
  energiedifferenz = e2 - e1;

  if(energiedifferenz < 0){
      spin = -spin; //akzeptiere Spinflip
  } else {
    //Vorschlagswahrscheinlichkeit
    uniform_real_distribution<double> distribution_2(0, 1);
    double p = distribution_2(generator);
    if(p < exp(-(1/kT)*energiedifferenz)){
      spin = -spin; //akzeptiere Spinflip
    }
  }
  return spin;
}

VectorXd randbedingungen(int position, MatrixXd r, double N){
    int n = sqrt(N);
    int spin_oben, spin_unten, spin_rechts, spin_links;
    //Randbedingungen
    //Werte sind Positionen in r, nicht die richtigen x- und y-Werte!
    if(position + n < N){ //Größe des Randes, Anzahl der Teilchen
      spin_rechts = r(2, position + n); //Muss in Summe immer gleich der Anzahl der Teilchen sein
    }
    else{
      spin_rechts = r(2, position - (N - n));
    }
    if(position - n >= 0){
      spin_links = r(2, position - n);
    }
    else{
      spin_links = r(2, position + (N - n));
    }
    if(position + 1 < N){
      spin_oben = r(2, position + 1);
    }
    else{
      spin_oben = r(2, position - n);
    }
    if(position - 1 >= 0){
      spin_unten = r(2, position - 1);
    }
    else{
      spin_unten = r(2, position + n);
    }
    VectorXd naechsteNachbarn = VectorXd::Zero(4);
    naechsteNachbarn << spin_oben, spin_unten, spin_rechts, spin_links;
    return naechsteNachbarn;
}

void aequilibrierung(MatrixXd &r, double S, VectorXd &energie, double N, double kT){
  VectorXd energie_spin, naechsteNachbarn;
  double Z;
  int t = 0, i;
  uniform_int_distribution<int> distribution(0, r.cols()-1); //nicht in while Schleife
  while(t < S){
    energie_spin = VectorXd::Zero(N);
    naechsteNachbarn = VectorXd::Zero(N);
    Z = 0;

    i = distribution(generator);
    r(2, i) = metropolis(r(2, i), randbedingungen(i, r, N), energie, i, kT);

    for(int k = 0; k < r.cols(); k++){
      naechsteNachbarn = randbedingungen(k, r, N);
      for(int j = 0; j < 4; j++){
        energie_spin(k) -= r(2, k)*naechsteNachbarn(j);
      }
    }

    for(int s = 0; s < r.cols(); s++){
      energie(t) += energie_spin(s);
    }
    energie(t) = energie(t)/N;
    t++;
  }
}

void spinflips(MatrixXd &r, double S, double N, double kT, VectorXd &magnetisierung, VectorXd &energie, VectorXd &betrag_magnetisierung){
  VectorXd energie_spin, magnetisierung_spin, betrag_magnetisierung_spin, naechsteNachbarn;
  double Z;
  int t = 0, i;
  uniform_int_distribution<int> distribution(0, r.cols()-1);
  while(t < S){ //Sweeps
    //Messungen nur nach jedem Sweep
    energie_spin = VectorXd::Zero(S);
    magnetisierung_spin = VectorXd::Zero(S);
    betrag_magnetisierung_spin = VectorXd::Zero(S);
    naechsteNachbarn = VectorXd::Zero(S);
    Z = 0;

    for(int j = 0; j < r.cols(); j++){ //Spinflips
      i = distribution(generator);
      r(2, i) = metropolis(r(2, i), randbedingungen(i, r, N), energie, i, kT);
    }

    for(int k = 0; k < r.cols(); k++){
      naechsteNachbarn = randbedingungen(k, r, N);
      for(int j = 0; j < 4; j++){
        energie_spin(k) -= r(2, k)*naechsteNachbarn(j);
      }
      magnetisierung_spin(k) += r(2, k);
      betrag_magnetisierung_spin(k) += r(2, k);
    }

    for(int s = 0; s < r.cols(); s++){
      energie(t) += energie_spin(s);
      magnetisierung(t) += magnetisierung_spin(s);
      betrag_magnetisierung(t) += betrag_magnetisierung_spin(s);
    }

    energie(t) = energie(t)/N;
    magnetisierung(t) = magnetisierung(t)/N;
    betrag_magnetisierung(t) = abs(betrag_magnetisierung(t))/N;

    t++;
  }
}

int main(){
  const int N = 10000;
  const int L = 100;
  double S = 5*1e4;
  int ordnung;

  double kT_1 = 1.5; //2.269185314
  double kT_2 = 2.269185314;
  double kT_3 = 3; //2.269185314

  MatrixXd r_1 = MatrixXd::Zero(3, N);
  MatrixXd r_2 = MatrixXd::Zero(3, N);
  MatrixXd r_3 = MatrixXd::Zero(3, N);

  VectorXd energie_1 = VectorXd::Zero(S);
  VectorXd energie_2 = VectorXd::Zero(S);
  VectorXd energie_3 = VectorXd::Zero(S);

  VectorXd magnetisierung_1 = VectorXd::Zero(S);
  VectorXd magnetisierung_2 = VectorXd::Zero(S);
  VectorXd magnetisierung_3 = VectorXd::Zero(S);

  VectorXd betrag_magnetisierung_1 = VectorXd::Zero(S);
  VectorXd betrag_magnetisierung_2 = VectorXd::Zero(S);
  VectorXd betrag_magnetisierung_3 = VectorXd::Zero(S);

  std::ofstream file;

  //Geordnete Anfangsbedingungen

  ordnung = 1;
  initialisierung(r_1, N, L, ordnung);
  aequilibrierung(r_1, S, energie_1, N, kT_1);

  file.open("Energie_Equi_geordnet_1.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < energie_1.size(); i++){
    file << energie_1(i) << "\n";
  }
  file.flush();
  file.close();

  initialisierung(r_2, N, L, ordnung);
  aequilibrierung(r_2, S, energie_2, N, kT_2);

  file.open("Energie_Equi_geordnet_2.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < energie_2.size(); i++){
    file << energie_2(i) << "\n";
  }
  file.flush();
  file.close();

  initialisierung(r_3, N, L, ordnung);
  aequilibrierung(r_3, S, energie_3, N, kT_3);

  file.open("Energie_Equi_geordnet_3.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < energie_3.size(); i++){
    file << energie_3(i) << "\n";
  }
  file.flush();
  file.close();

  //Zufällig angeordnete Spins

  ordnung = 0;

  r_1 = MatrixXd::Zero(3, N);
  r_2 = MatrixXd::Zero(3, N);
  r_3 = MatrixXd::Zero(3, N);

  energie_1 = VectorXd::Zero(S);
  energie_2 = VectorXd::Zero(S);
  energie_3 = VectorXd::Zero(S);

  initialisierung(r_1, N, L, ordnung);
  aequilibrierung(r_1, S, energie_1, N, kT_1);

  file.open("Energie_Equi_1.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < energie_1.size(); i++){
    file << energie_1(i) << "\n";
  }
  file.flush();
  file.close();

  spinflips(r_1, S, N, kT_1, magnetisierung_1, energie_1, betrag_magnetisierung_1);

  initialisierung(r_2, N, L, ordnung);
  aequilibrierung(r_2, S, energie_2, N, kT_2);

  file.open("Energie_Equi_2.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < energie_2.size(); i++){
    file << energie_2(i) << "\n";
  }
  file.flush();
  file.close();

  spinflips(r_2, S, N, kT_2, magnetisierung_2, energie_2, betrag_magnetisierung_2);

  initialisierung(r_3, N, L, ordnung);
  aequilibrierung(r_3, S, energie_3, N, kT_3);

  file.open("Energie_Equi_3.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < energie_3.size(); i++){
    file << energie_3(i) << "\n";
  }
  file.flush();
  file.close();

  spinflips(r_3, S, N, kT_3, magnetisierung_3, energie_3, betrag_magnetisierung_3);

  file.open("spinup_1.txt", std::ios_base::trunc);
  file << "#x \t y \n";
  for(int i = 0; i < r_1.cols(); i++){
    if(r_1(2, i) > 0){
      file << r_1(0, i) << "\t" << r_1(1, i) << "\n";
    }
  }
  file.flush();
  file.close();

  file.open("spindown_1.txt", std::ios_base::trunc);
  file << "#x \t y \n";
  for(int i = 0; i < r_1.cols(); i++){
    if(r_1(2, i) < 0){
      file << r_1(0, i) << "\t" << r_1(1, i) << "\n";
    }
  }
  file.flush();
  file.close();

  file.open("spinup_2.txt", std::ios_base::trunc);
  file << "#x \t y \n";
  for(int i = 0; i < r_2.cols(); i++){
    if(r_2(2, i) > 0){
      file << r_2(0, i) << "\t" << r_2(1, i) << "\n";
    }
  }
  file.flush();
  file.close();

  file.open("spindown_2.txt", std::ios_base::trunc);
  file << "#x \t y \n";
  for(int i = 0; i < r_2.cols(); i++){
    if(r_2(2, i) < 0){
      file << r_2(0, i) << "\t" << r_2(1, i) << "\n";
    }
  }
  file.flush();
  file.close();

  file.open("spinup_3.txt", std::ios_base::trunc);
  file << "#x \t y \n";
  for(int i = 0; i < r_3.cols(); i++){
    if(r_3(2, i) > 0){
      file << r_3(0, i) << "\t" << r_3(1, i) << "\n";
    }
  }
  file.flush();
  file.close();

  file.open("spindown_3.txt", std::ios_base::trunc);
  file << "#x \t y \n";
  for(int i = 0; i < r_3.cols(); i++){
    if(r_3(2, i) < 0){
      file << r_3(0, i) << "\t" << r_3(1, i) << "\n";
    }
  }
  file.flush();
  file.close();

  //file.open("Orte.txt", std::ios_base::trunc);
  //file << "#x \t y \n";
  //for(int i = 0; i < r.cols(); i++){
  //  file << r(0, i) << "\t" << r(1, i) << "\n";
  //}
  //file.flush();
  //file.close();

  file.open("Energie_1.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < energie_1.size(); i++){
    file << energie_1(i) << "\n";
  }
  file.flush();
  file.close();

  file.open("Energie_2.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < energie_2.size(); i++){
    file << energie_2(i) << "\n";
  }
  file.flush();
  file.close();

  file.open("Energie_3.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < energie_3.size(); i++){
    file << energie_3(i) << "\n";
  }
  file.flush();
  file.close();

  file.open("Magnetisierung_1.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < magnetisierung_1.size(); i++){
    file << magnetisierung_1(i) << "\n";
  }
  file.flush();
  file.close();

  file.open("Magnetisierung_2.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < magnetisierung_2.size(); i++){
    file << magnetisierung_2(i) << "\n";
  }
  file.flush();
  file.close();

  file.open("Magnetisierung_3.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < magnetisierung_3.size(); i++){
    file << magnetisierung_3(i) << "\n";
  }
  file.flush();
  file.close();

  file.open("Betrag_Magnetisierung_1.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < betrag_magnetisierung_1.size(); i++){
    file << betrag_magnetisierung_1(i) << "\n";
  }
  file.flush();
  file.close();

  file.open("Betrag_Magnetisierung_2.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < betrag_magnetisierung_2.size(); i++){
    file << betrag_magnetisierung_2(i) << "\n";
  }
  file.flush();
  file.close();

  file.open("Betrag_Magnetisierung_3.txt", std::ios_base::trunc);
  file << "#x \n";
  for(int i = 0; i < betrag_magnetisierung_3.size(); i++){
    file << betrag_magnetisierung_3(i) << "\n";
  }
  file.flush();
  file.close();

  return 0;
}