#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

double I_a(double x, double x_, double y_, double z_){
  return 1/sqrt((x-x_)*(x-x_) + y_*y_ + z_*z_);
}

double I_b(double x, double x_, double y_, double z_){
  return x_ / sqrt((x-x_)*(x-x_) + y_*y_ + z_*z_);
}

double mittelpunktsregel(double a, double b, double h, double x, double(*funptr)(double, double, double, double)){
  double integral = 0;
  double n = (b - a)/h;
  for(int i = 1; i < n; i++){
    for(int k = 1; k < n; k++){
      for(int j = 1; j < n; j++){
        integral += funptr(x, a + i*h - (h/2), a + k*h - (h/2), a + j*h - (h/2));
      }
    }
  }
  integral = h*h*h*integral;
  return integral;
}

//Intgerationsfunktion
//Übergeben werden die Intervallgrenzen, die verwendete Methode und die zu integrierende Funktion
double integration(double a, double b, double x, double(*funptr_2)(double, double, double, double)){
  double h = 1/(30*M_PI); //Start-Schrittweite
  double tmp_1; //Funktionswerte zwischenspeichern
  tmp_1 = mittelpunktsregel(a, b, h, x, funptr_2);
  return tmp_1; //Ergebnis zurückgeben
}

int main() {
  VectorXd x = 0.1*VectorXd::LinSpaced(70, 11, 80);

  // Aussen
  std::ofstream file;
  file.open("out.txt", std::ios_base::trunc);
  file << "# x phi_a(x) phi_b(x) \n";
  file << "x;phi_a;phi_b" << endl;
  for(int i = 0; i<x.size(); i++){
    file << x(i) << ";";
    file << integration(-1, 1, x(i), I_a) << ";";
    file << integration(-1, 1, x(i), I_b) << endl;
  }
  file.close();

  // Innen
  x = 0.1*VectorXd::LinSpaced(11, 0, 10);

  file.open("in.txt", std::ios_base::trunc);
  file << "# x phi_a(x) phi_b(x) \n";
  file << "x;phi_a;phi_b" << endl;
  for(int i = 0; i<x.size(); i++){
    file << x(i) << ";";
    file << integration(-1, 1, x(i), I_a) << ";";
    file << integration(-1, 1, x(i), I_b) << endl;
  }
  file.close();

  return 0;
}
