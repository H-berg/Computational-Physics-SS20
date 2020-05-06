#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "Eigen/Dense"
#include "./profiler.h"

using namespace std;
using namespace Eigen;

//Definition einer Funktion die Schritte 1-3 aus Aufgabenteil a ausführt und die Zeiten als Vektor zurückgibt
Eigen::VectorXd LUZerlegung(const int N)
{

  Vector3d t;
  Profiler::init(3);

  // Schritt 1: Erstellen einer zufälligen N x N Matrix sowie eines N-dimensionalen Vektors
  VectorXd b = VectorXd::Random(N);

  Profiler::start(0);
  MatrixXd M = MatrixXd::Random(N,N);
  Profiler::stop(0);
  t(0) = Profiler::getTimeInS(0);

  //Schritt 2: LU-Zerlegung 
  Profiler::start(1);
  PartialPivLU<MatrixXd> LU(M);
  Profiler::stop(1);
  t(1) = Profiler::getTimeInS(1);

  //Schritt 3: Lösen des Gleichungssystems
  Profiler::start(2);
  VectorXd x = LU.solve(b);
  Profiler::stop(2);
  t(2) = Profiler::getTimeInS(2);

  Profiler::resetAll();
  return t;

};



int main()
{
  const int n = 13; // Verwendung einer logarithmischen Skalierung; also ist NMax = 2^n
  int N_Max = pow(2, n);
  MatrixXd A(n+1, 3);
  int index = 0;

  for(int i = 1; i <= N_Max; i = 2*i)
  {
    VectorXd t = LUZerlegung(i);
    A.row(index) = t;
    index++;

  }

  ofstream f;
  f.open("data2b.txt");
  f << A << "\n";
  f.flush();
  f.close();

  return 0;
}
