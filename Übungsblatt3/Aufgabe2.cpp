#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

//using namespace std;
using namespace Eigen;
using namespace std;

//Die Funktion Federpendel gibt die Eigenwerte des Differentialgleichungssystems aus
VectorXd Federpendel(VectorXd &m, VectorXd &k, VectorXd &l, int N)
{
  //aufstellen der zu diagonalisierenden Matrix
  MatrixXd A(N, N);

  //die erste und letzte Zeile werden separat implementiert
  A(0, 0) = -k(0)/m(0);
  A(0, 1) = k(0)/m(0);
  A(N-1, N-1) = -k(N-2)/m(N-1);
  A(N-1, N-2) = k(N-2)/m(N-1);

  //Die Matrix wird mit den berechneten Elementen befüllt
  for(int i = 1; i <= N-2; i++)
  {
    A(i, i) = (-k(i)-k(i-1))/m(i);
    A(i, i-1) = k(i-1)/m(i);
    A(i, i+1) = k(i)/m(i);
  }

  //Die Eigenwerte der Matrix werden mittels eigen Bibliothek bestimmt und ausgegeben
  EigenSolver<MatrixXd> es(A);
  VectorXd eigenwerte = es.eigenvalues().real(); //Die Eigenwerte sind eh alle real, aber wir wollen nur eine Spalte haben
  return eigenwerte;

}


int main ()
{
  //Initialisierung der Anzahl der Massen N, sowie der Federkonstanten k, Massen m und Startwerte j
  const int N = 10;
  VectorXd m(N);
  VectorXd l(N-1);
  VectorXd k(N-1);

  for(int i = 0; i < N-1; i++)
  {
    m(i) = (i+1);
    k(i) = N-(i+1);
    l(i) = abs(5-(i+1))+1;
  }
  m(N-1) = N;

  //Berechnung der Eigenfrequenzen aus den Eigenwerten w = sqrt(-lambda)
  VectorXd eigenwerte = Federpendel(m, k, l, N);
  VectorXd eigenfrequenzen(N);

  for(int i = 0; i<N; i++)
  {
    eigenfrequenzen(i) = sqrt(abs(eigenwerte(i)));
  }

  ofstream f;
  f.open("Aufgabe2.txt");
  f << "Eigenfrequenzen" << "\n";
  f << eigenfrequenzen << "\n" ;
  f.flush();
  f.close();


  return 0;
}
