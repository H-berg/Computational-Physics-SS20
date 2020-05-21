#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

double I_1(double t){ //Integral (i)
    return exp(t)/t;
}

double H_I_1(double t){
  if(t == 0){
    return 2.0; //Grenzwert der Funktion ist 2
  } else{
    double r = (exp(t) - exp(-t))/t;
    return r;
  }
}


double I_2(double t){ //Integral (ii)
    return exp(-t)/sqrt(t);
}

double I_3(double t){ //Integral (ii)
  if(t == 0){
      return 1;
  } else{
      return sin(t)/t;
  }
}


//Übergeben werden jeweils die Intervallgrenzen, die Schrittweite und die zu integrierende Funktion
double trapezregel(double a, double b, double h, double(*funptr)(double)){
  double integral = 0;
  double sum = 0;
  double n = (b - a)/h;
  for(int i = 1; i < n ; i++){
    sum = sum + funptr(a + i*h);
  }
  integral = h*(funptr(a)/2 + funptr(b)/2) + h*sum;
  return integral;
}

double mittelpunktsregel(double a, double b, double h, double(*funptr)(double)){
  double integral = 0;
  double n = (b - a)/h;
  for(int i = 1; i < n ; i++){
    integral += funptr(a + i*h - (h/2));
  }
  integral = h*integral;
  return integral;
}

double simpsonregel(double a, double b, double h, double(*funptr)(double)){
  double integral = 0;
  double sum_1 = 0, sum_2 = 0;
  double n = (b - a)/h;
  for(int i = 1; i < n; i++){
    //double x = a + i * h;
    double k = funptr(a + i * h);
    if(i%2 == 0){
      sum_1 += k;
    }else{
      sum_2 += k;
    }
  }
  integral =  h/3*(funptr(a) + funptr(b) + 2*sum_1 + 4*sum_2);
  return  integral;
}

//Intgerationsfunktion
//Übergeben werden die Intervallgrenzen, die verwendete Methode und die zu integrierende Funktion
double integration(double a, double b, double(*funptr_1)(double, double, double, double(*funptr)(double)), double(*funptr_2)(double)){
  double h = 0.1; //Start-Schrittweite
  double tmp_1; //Funktionswerte zwischenspeichern
  double tmp_2;
  tmp_1 = funptr_1(a, b, h, funptr_2);
  tmp_2 = funptr_1(a, b, h/2, funptr_2);
  while(fabs(tmp_2 - tmp_1)/tmp_1 > 0.0005){ //Überprüfen der relativen Änderung
    h = h/2; //Schrittweite halbieren
    tmp_1 = funptr_1(a, b, h, funptr_2);
    tmp_2 = funptr_1(a, b, h/2, funptr_2);
  }
  return tmp_2; //Ergebnis zurückgeben
}

int main() {

  std::ofstream file;
  file.open("integration.txt", std::ios_base::trunc);
  file << "# \t trapezregel \n";
  //Das Integral wird in drei Teilintervalle aufgeteilt
  file << "I_1: \t" << integration(-1, -0.1, simpsonregel, I_1) + integration(0.1, 1, simpsonregel, I_1) + integration(-0.1, 0, simpsonregel, H_I_1)  << "\n";
  //Mittelpunktsregel da offene Newton-Cotes-Formel, Abbruch bei x_max = 100000
  file << "I_2: \t" << integration(0, 100000, mittelpunktsregel, I_2) << "\n";
  file << "I_3: \t" << 2*integration(0, 1000000, mittelpunktsregel, I_3) << "\n";

  file.close();

  return 0;
}
