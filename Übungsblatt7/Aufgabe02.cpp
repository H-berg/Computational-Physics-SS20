#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

VectorXd bewgl(VectorXd x, VectorXd x_p, double alpha){
  return -(x + alpha*x_p);
}

VectorXd runge_kutta(VectorXd (*func)(VectorXd,VectorXd, double), VectorXd y, VectorXd y_p, double alpha, double h){

int d = y.size();
VectorXd k1(d), k2(d), k3(d), k4(d);

k1 = h*func(y, y_p, alpha);
k2 = h*func(y+0.5*k1, y_p+0.5*k1, alpha);
k3 = h*func(y+0.5*k2, y_p+0.5*k2, alpha);
k4 = h*func(y+k3, y_p+k3, alpha);
return 1.0/6*(k1 + 2*k2 + 2*k3 + k4);
}


//  trägt die ersten d Einträge in den neuen Vektor ein, übergibt dann Werte
//  zu Berechnung der Beschleungigung. Dann wird Geschwindigkeit und Beschleungigung
//  zurück gegeben.
VectorXd swap(VectorXd y, double alpha){
  int d = y.size();
  VectorXd a(d);
  d *= 0.5;

  a.segment(0,d) = y.segment(d,d);
  a.segment(d,d) = bewgl(y.segment(0,d), y.segment(d,d),alpha);
  return a;
}

void adams_bashfort(VectorXd (*func)(VectorXd, double),
                                  double alpha,
                                  MatrixXd &Y,
                                  double &h,
                                  int n){

  Y.col(n+1) = 55*func(Y.col(n), alpha);
  Y.col(n+1) -= 59*func(Y.col(n-1), alpha);
  Y.col(n+1) += 37*func(Y.col(n-2), alpha);
  Y.col(n+1) -= 9*func(Y.col(n-3), alpha);
  Y.col(n+1) *= h/24.0;
  Y.col(n+1) += Y.col(n);
}


int main(){
  const int N = 300, T = 20, d=3;
  double alpha = 0.1, h = T/double(N);
  MatrixXd kutta(2*d, N+2);
  VectorXd r(d), v(d), Energie(N+2);
  r << 1, 2, 3;
  v << 0, 1, 0;

// Aufgabe a) und b) ----------------------------------------------------------

  //4 Startwerte mittels Kutta für alpha = 0.1
  for(int i=0; i < 4; i++){
    kutta.col(i).segment(0,d) = r.transpose();
    kutta.col(i).segment(d,d) = v.transpose();

    Energie(i) = kutta.col(i).segment(0,d).dot(kutta.col(i).segment(0,d));
    Energie(i) += kutta.col(i).segment(d,d).dot(kutta.col(i).segment(d,d));
    Energie(i) *= 0.5;

    v = v + runge_kutta(bewgl, r, v, alpha, h);
    r = r + runge_kutta(bewgl, r, v, 1, h);
  }

  for(int i = 3; i <= N; i++){
    adams_bashfort(swap, alpha, kutta, h, i);
    Energie(i+1) = kutta.col(i).segment(0,d).dot(kutta.col(i).segment(0,d));
    Energie(i+1) += kutta.col(i).segment(d,d).dot(kutta.col(i).segment(d,d));
    Energie(i+1) *= 0.5;
  }

  // Die Geschwindigkeiten brauchen wir aber hier nicht mehr
  // Zeit speichern in kutta Matrix.
  kutta.row(3) = VectorXd::LinSpaced(N+2, 0, T);

  // Energie ersetzt den Platz in Kutta für die Geschwindigkeit.
  kutta.row(4) = Energie;


  std::ofstream file;
  file.open("A2_1a.txt", std::ios_base::trunc);
  file << "# x y z t E \n";
  file << kutta.transpose().leftCols(5) << endl;
  file.close();

  alpha = 0;

  //4 Startwerte mittels Kutta
  for(int i=0; i < 4; i++){
    kutta.col(i).segment(0,3) = r.transpose();
    kutta.col(i).segment(d,d) = v.transpose();

    v = v + runge_kutta(bewgl, r, v, alpha, h);
    r = r + runge_kutta(bewgl, r, v, 1, h);
  }

  for(int i = 3; i <= N; i++){
    adams_bashfort(swap, alpha, kutta, h, i);
  }


  file.open("A2_2a.txt", std::ios_base::trunc);
  file << "# x y z \n";
  file << kutta.transpose().leftCols(3) << endl;
  file.close();

  alpha = -0.1;

  //4 Startwerte mittels Kutta
  for(int i=0; i < 4; i++){
    kutta.col(i).segment(0,3) = r.transpose();
    kutta.col(i).segment(d,d) = v.transpose();

    v = v + runge_kutta(bewgl, r, v, alpha, h);
    r = r + runge_kutta(bewgl, r, v, 1, h);
  }

  for(int i = 3; i <= N; i++){
    adams_bashfort(swap, alpha, kutta, h, i);
  }


  file.open("A2_3a.txt", std::ios_base::trunc);
  file << "# x y z \n";
  file << kutta.transpose().leftCols(3) << endl;
  file.close();


}
