#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include "profiler.h"
using namespace std;
using namespace Eigen;
using Eigen::PartialPivLU;


int main() {
  //Aufgabenteil a
  Profiler::init(6);

  double inv;
  double fullpiv;
  double partialpiv;

  std::ofstream file;
  file.open("times.txt", std::ios_base::trunc);
  file << "#n \t inv \t fullpiv \t partialpiv \t A_to_B \t A_to_C \t B_to_C \n";

  for(int i=1; i<=1000; i++){
    Profiler::resetAll();

    MatrixXd M = MatrixXd::Random(i, i);
    VectorXd b = VectorXd::Random(i);

    //Lösung mit inverser Matrix
    Profiler::start(0);
    MatrixXd Inv = M.inverse();
    VectorXd A = Inv*b;
    Profiler::stop(0);

    //Lösung mit voller Pivotisierung
    Profiler::start(1);
    VectorXd B = M.fullPivLu().solve(b);
    Profiler::stop(1);

    //Lösung mit teilweiser Pivotisierung
    Profiler::start(2);
    VectorXd C = M.partialPivLu().solve(b);
    Profiler::stop(2);

    //Werte sollen nicht gespeichert werden, falls sie ungleich sind
    if(A.isApprox(B, 1E-8) != true || A.isApprox(C, 1E-8) != true || B.isApprox(C, 1E-8) != true || M.determinant() == 0){
      continue;
    }

    inv = Profiler::getTimeInS(0);
    fullpiv = Profiler::getTimeInS(1);
    partialpiv = Profiler::getTimeInS(2);

    file << i << "\t" << inv << "\t" << fullpiv << "\t" << partialpiv << "\n";

  }

  file.flush();
  file.close();

  return 0;
}
