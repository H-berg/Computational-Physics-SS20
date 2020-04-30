#include <iostream>
#include <fstream>
#include <eigen3/Eigen/LU>

using namespace std;
using namespace Eigen;

int main(){
  //b)------------------------------------------------------------------------
  //cout << "b)" << endl;

  //cout << "Das zu lösende lin. GS: Ax'=x" << std::endl;

  // Vektor x=(2,0,2)
  Vector3d x;
  x << 2,0,2;
  //cout << "Vektor x:\n" << x << endl;

  // Matrix A aus den Basen a_1, a_2, a_3
  Matrix3d A;
  A << 0.5,-0.5,0,
       sqrt(3)/2,sqrt(3)/2,0,
       0,0,1;
  //cout << "Matrix A:\n" << A << endl;

  // da sich A durch PartialPivLU ändert, wird eine kopie erstellt
  Matrix3d A0 = A;

  // LU-Zerlegung von Matrix A durch Eigen:PartialPivLu
  //PartialPivLU<MAtrix3d> dec(A);
  PartialPivLU<Ref<Matrix3d> > lu(A0);

  // Permutations-Matrix P
  Matrix3d P = lu.permutationP();

  //cout << "Permutations-Matrix P:\n" << P << endl;
  //LU-Martix, welche aus Lower- und Upper-triangular-Matrix besteht
  Matrix3d L = lu.matrixLU().triangularView<Lower>();
  //Für größere Matrizen mit Schleife, hier sollte es so gehen
  L(0,0)=1;
  L(1,1)=1;
  L(2,2)=1;

  //U-Matrix bestimmen
  Matrix3d U = lu.matrixLU().triangularView<Upper>();

  //Vektor x' bestimmen
  Vector3d x_strich = lu.solve(x);
  // teste zunächst ob x' richtig ist. Danach Ausgabe von x'
  //cout << "Test: A*x'-x=" << (A * x_strich - x).norm() << endl;
  //cout << "LU-Matrix :\n" << x_strich << endl;

  //einmal abspeichern
  fstream f;
  f.open("data1b.txt", ios::out);
  f << "P=\n"<< P << "\n\n";
  f << "L=\n"<< L << "\n\n";
  f << "U=\n"<< U << "\n\n";
  f << "x'=\n"<< x_strich << "\n\n";
  f << "Test: A*x'-x=" << (A * x_strich - x).norm() << "\n\n";

  f.flush();
  f.close();

  //c)------------------------------------------------------------------------
  //cout << "c)" << endl;

  //cout << "Matrix A:\n" << A << endl;

  Vector3d y;
  y <<1,2*sqrt(3),3;
  Vector3d y_strich = lu.solve(y);
  // teste zunächst ob y' richtig ist. Danach Ausgabe von y'
  //cout << "Test: A*x'-x=" << (A * y_strich - y).norm() << endl;
  //cout << "LU-Matrix :\n" << y_strich << endl;

  // abspeichern
  f.open("data1c.txt", ios::out);
  f << "y'=\n"<< y_strich << "\n\n";

  f.flush();
  f.close();
  //d)------------------------------------------------------------------------
  //cout << "d)" << endl;

  Matrix3d B;
  B << 0,-0.5,0.5,
       0,sqrt(3)/2,sqrt(3)/2,
       1,0,0;

  //cout << "Matrix B:\n" << B << endl;
  // erneuerte LU-Zerlegung von Matrix B durch Eigen:PartialPivLu
  lu.compute(B);

  // Permutations-Matrix P
  Matrix3d P_B = lu.permutationP();
  //cout << "Permutations-Matrix P:\n" << P_B << endl;
  //LU-Martix, welche aus Lower- und Upper-triangular-Matrix besteht

  Matrix3d L_B = lu.matrixLU().triangularView<Lower>();
  //Für größere Matrizen mit Schleife, hier sollte es so gehen
  L_B(0,0)=1;
  L_B(1,1)=1;
  L_B(2,2)=1;

  //U-Matrix bestimmen
  Matrix3d U_B = lu.matrixLU().triangularView<Upper>();

  //abspeichern
  f.open("data1d.txt", ios::out);
  f << "P_B=\n"<< P_B << "\n\n";
  f << "L_B=\n"<< L_B << "\n\n";
  f << "U_B=\n"<< U_B << "\n\n";

  f.flush();
  f.close();
}
