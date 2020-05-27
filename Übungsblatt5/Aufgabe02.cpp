#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <tuple>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

double func(double x){
  return x*x-2;
}

// 1. Ableitung mit O(h^2) Skritp4 S.1 4.1A)
double diff1(double (*funptr)(double), double x){
  double h = 1e-3;
  return 0.5*(funptr(x+h)-funptr(x-h))/h;
}
// 2. Ableitung mit O(h^2)  Skritp4 S.1 4.1B)
double diff2(double (*funptr)(double), double x){
  double h = 1e-2;
  return (funptr(x+h) - 2*funptr(x) + funptr(x-h))/(h*h);
}

// Intervallhalbierung nach Kirfeld-Skript
int ih(double (*funptr)(double), double a, double b, double c, double gs){
  double u;
  int i = 0;
  while((c-a) > gs)
  {
    // halbiere das längere Intervall
    if(abs(c-b) > abs(b-a))
    {
      u = (b+c)*0.5;
      if(funptr(u) < funptr(b))
      {
        a = b;
        b = u;
      }
      else
      {
        c = u;
      }
    }

    else
    {
      u = (a+b)*0.5;
      if(funptr(u) < funptr(b))
      {
        c = b;
        b = u;
      }
      else
      {
        a = u;
      }
    }
    i++;
  }
  return i;
}

// Newton nach Skript 5.1.3
int newton(double (*funptr)(double), double x0, double gs){
  double h;
  int i = 0;
  do
  {
    h = diff1(funptr, x0)/diff2(funptr, x0);
    x0 = x0 - h;
    i++;
  }
  while(abs(h) > gs);

  return i;
}

int main(){
  std::ofstream file;
  file.open("A2.txt", std::ios_base::trunc);
  file << "Die Intervallhalbierung benötigt " << ih(func,-0.5,-0.1, 2,1e-9)
  << " Interrationsschritte." << endl;
  file << "Das Newton-Verfahren benötigt " << newton(func,1,1e-9)
  << " Interrationsschritte." << endl;
  file.close();
}
