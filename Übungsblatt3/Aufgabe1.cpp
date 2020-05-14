#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

//using namespace std;
using namespace Eigen;
using namespace std;
const int N = 50;


//Funktion zur Berechnung des kleinsten Eigenwerts
double lanczos(MatrixXd &A, VectorXd q){

  //Initialisierung
  MatrixXd one = MatrixXd::Identity(N, N);
  VectorXd q0 = VectorXd::Zero(N);
  VectorXd qtemp;
  VectorXd v(N);
  MatrixXd T(1, 1);
  VectorXd eigenwerte;
  EigenSolver<MatrixXd> es;
  VectorXd::Index position;
  int i = 0;


  double delta = (q.transpose()*A)*q;
  double gamma, gamma0, gammatemp;
  double EW_then = 1000;
  double EW_now = 0;
  gamma0 = 1;
  T(0, 0) = delta;



  while(gamma != 0 && abs(EW_then-EW_now) > 0.001){
    //Iteration
    v = (A-delta*one)*q - gamma0*q0;

    gammatemp = gamma;
    gamma = sqrt(v.dot(v));
    gamma0 = gammatemp;

    qtemp = q;
    q = 1/gamma*v;
    q0 = qtemp;

    delta = (q.transpose()*A)*q;

    T.conservativeResize(i+2, i+2); // die Größe der Matrix wird verändert ohne die Elemente zu beeinflussen
    T(i+1, i) = gamma;
    T(i, i+1) = gamma;
    T(i+1, i+1) = delta;
    i++;

    //Untersuchung der Eigwenwerte
    EW_then = EW_now;
    es.compute(T, true);
    eigenwerte = es.eigenvalues().real();

    //Wir wollen das Betragmäßig kleinste Element, dafür gibt es bestimmt einen Befehl aber ich mach das mal so
    for(int j = 0; j < eigenwerte.size(); j++)
    {
      eigenwerte(j) = abs(eigenwerte(j));
    }

    //Finde den kleinsten Eigenwert und schreibe dessen Index in position
    EW_now = eigenwerte.minCoeff(&position);
  }

  //Abspeichern des Eigenwertes und Eigenvektors in einen Vektor
  return EW_now;
}


int main()
{
  const double t = 1;
  const double epsilon = -5;

  MatrixXd one = MatrixXd::Identity(N, N);


  //Erstellen eines zufälligen normierten Vektors
  VectorXd q = VectorXd::Random(N);
  q = q/(sqrt(q.dot(q)));

  MatrixXd H = MatrixXd::Zero(N, N);

  for(int i = 0; i<N-1; i++)
  {
    H(i+1, i) = -t;
    H(i, i+1) = -t;
  }

  H(0, N-1) = -t;
  H(N-1, 0) = -t;
  H(N/2, N/2) = epsilon;
  double eigenwert = lanczos (H, q);

  MatrixXd M = H - one*eigenwert;
  //Eigenvektor berechnen
  FullPivLU<MatrixXd> LUf(M);
  VectorXd eigenvektor = LUf.solve(VectorXd::Zero(N));

  ofstream f;
  f.open("Aufgabe1b_EV_EW.txt");
  f << "Der Eigenwert lautet: " << eigenwert << "\n";
  f << "Der Eigenvektor lautet" << eigenvektor << "\n";
  f.flush();
  f.close();



  return 0;


}
