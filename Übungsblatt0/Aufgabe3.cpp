#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// Schrittweite für beide Verfahren, global einstellbar
double deltat = 0.01;

// gibt den Wert von y an der Stelle t=x aus für y'=-y
// Euler-Verfahren
double euler(double x, double y0)  // Startparameter werden mitübergeben, da b) so einfacher lösbar ist
{
  //double y0 = 1;
  double y=1;
  // double deltat = 0.1;
  int i = 0;
  while(i*deltat<x)
  {
    y=y0*(1-deltat);
    y0 = y;
    i++;
  }
  return y;
}

// symmetrischen Euler-Verfahren
double symEuler(double x, double y0, double y1) // Startparameter werden mitübergeben, da b) so einfacher lösbar ist
{
  //double y0 = 1;
  // double deltat = 0.1;
  //double y1 = exp(-deltat);
  double y = 1;
  int i = 0;
  while(i*deltat<x)
  {
    y=-2*deltat*y1+y0;
    y0 = y1;
    y1 = y;
    i++;
  }
  return y;
}

// analytische Funktion y(t)=exp(-t)
double analytisch(double x)
{
  return exp(-x);
}

int main()
{
  //x-Wert
  double x;
  //Setzen der Starparameter für a)
  double y0=1;
  double y1 = exp(-deltat);

  fstream f;
  f.open("data3a.txt", ios::out);
  for(int i=0; i<100; i++){
    x = i*0.1; // x nimmt Werte zwischen 0 und 10 in 0.1-Schritten an
    f << x << "\t\t\t" <<euler(x, y0) << "\t\t\t" << symEuler(x, y0, y1)<< "\t\t\t" << analytisch(x)<<"\n"; // Frage: Wie kann ich die Ausgabe 'schön' untereinander schreiben?
  }
  f.flush();
  f.close();

  //Setzen der Starparameter für b) ansonsten ist das Verfahren das gleiche wie in a)
  y0=1-deltat;
  y1= 1-deltat;
  f.open("data3b.txt", ios::out);
  for(int i=0; i<100; i++){
    x = i*0.1; // x nimmt Werte zwischen 0 und 10 in 0.1-Schritten an
    f << x << "\t\t\t" <<euler(x, y0) << "\t\t\t" << symEuler(x, 1, y1)<< "\t\t\t" << analytisch(x)<<"\n";
  }
  f.flush();
  f.close();

  return 0;
}
