#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>

#include "moodle/service.cpp"

using namespace std;
using namespace Eigen;

// zu b) Approximationsfunktion -----------------------------------------------
MatrixXd Approximation(MatrixXd &mat, int k, MatrixXd &U, MatrixXd &V, VectorXd &W){
  MatrixXd M(mat.rows(), mat.cols());
  MatrixXd temp(V.cols(), V.rows());
  temp = V.transpose();
  
  if(k==0){
    M = mat;
  }
  else{
    for (int i = 0; i < k; i++)
    {
      M += W(i) * U.col(i) * temp.row(i);
    }
  }
  return M;
};

int main(){
  // a) Daten einlesen und SVD durchführen ------------------------------------
  MatrixXd A, A2;
  loadData(A2, "moodle/Bild", 512, 512);

  A = A2.transpose(); // Die Matrix wird transponiert, da das Bild sonst falsch herum ist

  BDCSVD<MatrixXd> svd(A, ComputeFullU|ComputeFullV);
  VectorXd sing = svd.singularValues();
  MatrixXd U = svd.matrixU();
  MatrixXd V = svd.matrixV();


  ////b) Rang-k-Approximation ---------------------------------------------------
  int k[] = {0,10,20,50};
  ofstream f;

  for(int i = 0; i < sizeof(k)/sizeof(k[0]); i++){
    f.open("bild_rang_"+to_string(k[i])+".txt");
    f << Approximation(A, k[i], U, V, sing);
    f.flush();
    f.close();
  }

  return 0;
}
