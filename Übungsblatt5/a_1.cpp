#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include <complex>

using namespace std;
using namespace Eigen;

//Gespiegelte Binärdarstellung der einzelnen Zahlen
int bin(int n, int l)
{
  VectorXd bin = VectorXd::Zero(l);
  int p = 0;
  do
  {
    if (n % 2 == 0)
    {
      bin(p) = 0;
      p++;
      n = n / 2;
    }
    else if (n % 2 == 1)
    {
      if (n > 1)
      {
        bin(p) = 1;
        p++;
        n = (n - 1) / 2;
      }
      else
      {
        n = 0;
        bin(p) = 1;
      }
    }
  } while (n > 0);
  int k = 0;
  for (int i = 0; i < bin.size(); i++)
  {
    if (bin(i) == 0)
    {
      k = k * 2;
    }
    else
    {
      k = k * 2 + 1;
    }
  }
  return k;
}

//Diskrete Fouriertransformation
VectorXcd DFT(int m, int n, VectorXcd f)
{
  const complex<double> I(0.0, 1.0);
  VectorXcd F = VectorXcd::Zero(n);
  MatrixXcd S_1 = MatrixXcd::Zero(n, n);
  MatrixXcd S_0 = MatrixXcd::Zero(n, n);
  VectorXcd f_l = VectorXcd::Zero(n);
  //Schritt 1
  for (int i = 0; i < n; i++)
  {
    int k = bin(i, m);
    f_l(i) = f(k);
  }
  for (int j = 0; j < n; j++)
  {
    for (int l = 0; l < n; l++)
      S_0(j, l) = f_l(l);
  }
  //Schritt 2
  for (int k = 1; k <= m; k++)
  {
    for (int l = 0; l <= pow(2, m - k) - 1; l++)
    {
      for (int j = 0; j <= pow(2, k) - 1; j++)
      {
        S_1(j, l) = S_0(j, 2 * l) + (cos(2 * (M_PI * j) / pow(2, k)) + I * sin(2 * (M_PI * j) / pow(2, k))) * S_0(j, 2 * l + 1);
        for (int i = 1; i <= pow(2, m - k) - 1; i++)
        {
          S_1(j + i * pow(2, k), l) = S_1(j, l);
        }
      }
    }
    S_0 = S_1;
  }
  for (int j = 0; j < n; j++)
  {
    F(j) = S_1(j, 0);
  }
  return F;
}

//Funktion für Aufgabenteil 1
VectorXcd f(int m, int n)
{
  VectorXcd f = VectorXcd::Zero(n);
  for (int l = 0; l < n; l++)
  {
    f(l) = sqrt(1 + l);
  }
  return f;
}

//Funktion für Aufgabenteil 2
double f_1(double x)
{
  return exp(-(x * x) / 2);
}

//Funktion für Aufgabenteil 3 (Rechteckschwingung)
double f_2(double x)
{
  if (x > 0 && x < M_PI)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}

//Schnelle Fouriertransformation mit Verschiebung und Phasenfaktormultiplikation
VectorXcd FFT(int m, int n, VectorXd x, double (*funptr)(double)){
  VectorXd H = VectorXd::Zero(n);
  VectorXcd F = VectorXcd::Zero(n);

  for (int i = 0; i < x.size(); i++){
    H(i) = funptr(x(i));
  }

  F = DFT(m, n, H);

  //Verschiebung der Elemente
  VectorXcd G = F;
  for (int i = 0; i < n / 2; i++){
    F(i) = G(i + n / 2);
    F(i + n / 2) = G(i);
  }

  //Multiplikation mit Phasenfaktor
  const complex<double> I(0.0, 1.0);
  double L = x(0) - x(x.size() - 1);
  double d_x = L / n;
  double x_min = x(0);

  for (int j = 0; j < n; j++){
    F(j) = F(j) * (d_x / (2 * M_PI)) * (cos((2 * M_PI * x_min * j) / L) + I * sin((2 * M_PI * x_min * j) / L));
  }
  return F;
}

//Direkte Berechnung
VectorXcd D(int m, int n, VectorXcd f){
  const complex<double> I(0.0, 1.0);
  VectorXcd F = VectorXcd::Zero(n);
  for (int j = 0; j < n; j++){
    for (int l = 0; l < n; l++){
      F(j) += (cos(2 * M_PI * j * l / n) + I * sin(2 * M_PI * j * l / n)) * sqrt(1 + l);
    }
  }
  return F;
}

int main(){
  int m = 7;
  int n = pow(2, m);
  VectorXd x = VectorXd::LinSpaced(n, -10, 10);     //x-Werte von Aufgabenteil 2
  VectorXd y = VectorXd::LinSpaced(n, -M_PI, M_PI); //x-Werte von Aufgabenteil 3

  //Aufgabenteil 1
  std::ofstream file;
  file.open("dft_3.txt", std::ios_base::trunc);
  file << "#DFT \n";
  file << DFT(3, 8, f(3, 8)) << "\n";
  file.flush();
  file.close();
  file.open("direkt_3.txt", std::ios_base::trunc);
  file << "#direkt \n";
  file << D(3, 8, f(3, 8)) << "\n";
  file.flush();
  file.close();
  file.open("dft_4.txt", std::ios_base::trunc);
  file << "#DFT \n";
  file << DFT(4, 16, f(4, 16)) << "\n";
  file.flush();
  file.close();
  file.open("direkt_4.txt", std::ios_base::trunc);
  file << "#direkt \n";
  file << D(4, 16, f(4, 16)) << "\n";
  file.flush();
  file.close();

  //Aufgabenteil 2 und 3
  file.open("fft_1.txt", std::ios_base::trunc);
  file << "#fft_1 \n";

  VectorXcd fft1 = FFT(m, n, x, f_1);
  VectorXcd fft2 = FFT(m, n, y, f_2);
  for (int i = 0; i < n; i++)
  {
    file << -64/(2*M_PI*(2*M_PI/x.size()))+i/(2*M_PI*(2*M_PI/x.size())) << "\t" << abs(fft1(i))<< "\n";
  }
  file.flush();
  file.close();
  file.open("fft_2.txt", std::ios_base::trunc);
  file << "#fft_2 \n";
  for (int i = 0; i < n; i++)
  {
    //Einlesen des Betrags
    file << abs(fft2(i)) << "\n";
  }
  file.flush();
  file.close();
}
