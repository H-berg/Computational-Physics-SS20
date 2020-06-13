#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <tuple>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
// Kiefeld Skript S. 55

// -------------------------------- 3-D --------------------------------------


// Kraftfeld F(r) = -r*m*w^2 (hier ist w = 1)
Vector3d kraftfeld(Vector3d r, double m){
  return -m*r;
}


Vector3d geschw(Vector3d r, double m){
  return r;
}


Vector3d runge_kutta(Vector3d (*func)(Vector3d, double), Vector3d y, double m, double &h){

// Zeitschritte nach Kierfeld-Skript S.50
const double M = 1/m;
Vector3d k1, k2, k3, k4;

k1 = h*M*func(y, m);
k2 = h*M*func(y+0.5*k1, m);
k3 = h*M*func(y+0.5*k2, m);
k4 = h*M*func(y+k3, m);
return 1.0/6*(k1 + 2*k2 + 2*k3 + k4);
}


int main(){
  const double m = 1;
  int N = 300;
  double T=20, h = T/N;

  Vector3d r, v;
  MatrixXd res(N+1,5);

//  Aufgabe 1a) --------------------------------------------------------------

  v << 0,0,0;
  r << 1, 2, 3;


  for(int i=0; i <= N; i++){
    res.row(i).segment(0,3) = r.transpose();
    res(i,3) = i*h;
    res(i,4) = 0.5*m*(v.dot(v) + r.dot(r));
    v = v + runge_kutta(kraftfeld, r, m, h);
    r = r + runge_kutta(geschw, v, 1, h);
  }

  std::ofstream file;
  file.open("teil_a1.txt", std::ios_base::trunc);
  file << "# v_dot = -r --> v kuttan --> r kuttan \n";
  file << res << endl;
  file.close();

  v << 1, 1, 1;
  r << 1, 2, 3;

  for(int i=0; i <= N; i++){
    res.row(i).segment(0,3) = r.transpose();
    res(i,3) = i*h;
    res(i,4) = 0.5*m*(v.dot(v) + r.dot(r));
    v = v + runge_kutta(kraftfeld, r, m, h);
    r = r + runge_kutta(geschw, v, 1, h);
  }

  file.open("teil_a2.txt", std::ios_base::trunc);
  file << "# v_dot = -r --> v kuttan --> r kuttan \n";
  file << res << endl;
  file.close();

  //  Aufgabe 1b) -------------------------------------------------------------
  N = 10;
  v << 1, 0, 3;
  r << 0, 1, 0;
  Vector3d ri = r;
  int p=0;

    h = 1e-8;     // durch händisches ausprobieren. In einer While-Schleife hat
                  // das irgendwie immer den gleichen Wert ausgespuckt :/
    for(int i=0; i < N; i++){
      v = v + runge_kutta(kraftfeld, ri, m, h);
      ri = ri + runge_kutta(geschw, v, 1, h);
    }


//  Aufgabe 1c) ---------------------------------------------------------------
r << 1, 2, 3;
v << 0, 0, 0;

T = 20;
N = 500;
h = T/N;

file.open("teil_c.txt", std::ios_base::trunc);
file << "# v_dot = -r --> v kuttan --> r kuttan \n";
for(int i=0; i <= N; i++){
  file << 0.5*m*(v.dot(v) + r.dot(r)) << "\t";
  file << i*h << endl;
  v = v + runge_kutta(kraftfeld, r, m, h);
  r = r + runge_kutta(geschw, v, 1, h);
}
file.close();


return 0;
}
