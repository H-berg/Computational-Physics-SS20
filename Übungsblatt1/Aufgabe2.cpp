#include <iostream>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"

//using namespace std;
//Frage: Kann es zu Problemen kommen, wenn man zwei using namespace verwendet?
using namespace Eigen;

int main(){

  //Vektoren erstellen.
  VectorXf x(10);
  VectorXf y(10);
  x << 0, 2.5, -6.3, 4, -3.2, 5.3, 10.1, 9.5, -5.4, 12.7;
  y << 4, 4.3, -3.9, 6.5, 0.7, 8.6, 13, 9.9, -3.6, 15.1;

  //Matrix erstellen 
  MatrixXf A(10, 2);
  A.col(0) = x;
  A.col(1) << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

  //Überführung in ein quadratisches Problem (Aufgabenteil b)
  MatrixXf P(2, 2);
  MatrixXf AT(2, 10);
  AT = A.transpose();
  P = AT*A;

  // das neue quadratische Problem lautet Pn=y2 (mit n =(m, n))
  Vector2f y2 = AT*y;

  //Lösen mittels LU-Zerlegung
  Vector2f n = P.partialPivLu().solve(y2);

  //abspeichern
  std::ofstream f;
  f.open("data2c.txt");
  f << "# m n\n" << n(0) << "\t" << n(1) << "\n";
  f.flush();
  f.close();

  //für den Plot werden auch noch einmal x- und y-Werte abgespeichert
  f.open("data2a.txt");
  f << "# x y \n";
  for(int i = 0; i<10; i++){
    f << x(i) << " \t" << y(i) << "\n";
  }

  f.flush();
  f.close();

  return 0;
}
