#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

double rosenbrock(double x_1, double x_2){
  return (1-x_1)*(1-x_1) + 100*(x_2-x_1*x_1)*(x_2-x_1*x_1);
}

//Definition der Ableitungen für die Komponenten des Gradienten
double ableitung_1(double x_1, double x_2, double(*f)(double, double)){
  double h = 1e-6;
  return (f(x_1+h, x_2)-f(x_1-h, x_2))/(2*h);
}

double ableitung_2(double x_1, double x_2, double(*f)(double, double)){
  double h = 1e-6;
  return (f(x_1, x_2+h)-f(x_1, x_2-h))/(2*h);
}

//Berechnung des Gradienten
VectorXd gradient(VectorXd x, double(*f)(double, double)){
  VectorXd g(x.size());
  g(0) = ableitung_1(x(0), x(1), f);
  g(1) = ableitung_2(x(0), x(1), f);
  return g;
}

//Newton-Verfahren (Verfahren zur Nullstellenbestimmung: Suche nach Nullstelle der ersten Ableitung)
double newton(VectorXd x, VectorXd g, double lambda, double(*f)(double, double)){
  double h = 1e-6;
  double d_x = 0, diff_1 = 0, diff_2 = 0;
  double lambda_2 = 0;
  do{
    diff_1 = (f(x(0)+g(0)*(lambda_2+h), x(1)+g(1)*(lambda_2+h))-f(x(0)*g(0)*(lambda_2-h), x(1)*g(1)*(lambda_2-h)))/(2*h);
    diff_2 = (f(x(0)+g(0)*(lambda_2+h), x(1)+g(1)*(lambda_2+h))-2*f(x(0)+g(0)*(lambda_2), x(1)+g(1)*(lambda_2))+f(x(0)+g(0)*(lambda_2-h), x(1)+g(0)*(lambda_2-h)))/(h*h);
    d_x = diff_1/diff_2;
    lambda_2 = lambda - d_x;
    lambda = lambda_2;
  } while(abs(d_x) > 1e-5);
  return lambda_2;
}

//Bestimmung der initialen Hesse-Matrix
//Verfahren 1: Exakte inverse Hesse-Matrix
MatrixXd H_exakt(VectorXd x0, double(*f)(double, double)){
  MatrixXd H = MatrixXd::Zero(2, 2);
  double tmp1 = 0;
  double tmp2 = 0;
  double h = 1e-6;

  H(0, 0) = (f(x0(0)+h, x0(1))-2*f(x0(0), x0(1))+f(x0(0)-h, x0(1)))/(h*h);
  tmp1 = ableitung_1(x0(0), x0(1), f);
  H(0, 1) = ableitung_2(tmp1, x0(1), f);
  tmp2 = ableitung_2(x0(0), x0(1), f);
  H(1, 0) = ableitung_1(tmp2, x0(1), f);
  H(1, 1) = (f(x0(0), x0(1)+h)-2*f(x0(0), x0(1))+f(x0(0), x0(1)-h))/(h*h);

  H = H.inverse();
  //cout << H << endl;
  return H;
}

//Verfahren 2: Diagonalmatrix
MatrixXd H_diagonal(VectorXd x0, double(*f)(double, double)){
  MatrixXd H = MatrixXd::Zero(2, 2);
  MatrixXd C = MatrixXd::Zero(2, 2);
  double h = 1e-6;

  H(0, 0) = (f(x0(0)+h, x0(1))-2*f(x0(0), x0(1))+f(x0(0)-h, x0(1)))/(h*h);
  H(1, 1) = (f(x0(0), x0(1)+h)-2*f(x0(0), x0(1))+f(x0(0), x0(1)-h))/(h*h);
  H.inverse();
  C(0, 0) = H(0, 0);
  C(1, 1) = H(1, 1);
  //cout << C << endl;
  return C;
}

//Verfahren 3: Vielfaches der Einheitsmatrix
MatrixXd H_einheitsmatrix(VectorXd x0, double(*f)(double, double)){
  MatrixXd I = MatrixXd::Identity(2, 2);
  MatrixXd H = MatrixXd::Zero(2, 2);
  double w;

  w = f(x0(0), x0(1));
  H = w * I;
  //cout << H << endl;
  return H;
}

VectorXd BFGS(VectorXd x0, MatrixXd c0, double epsilon, double(*f)(double, double), VectorXd(*gradient)(VectorXd, double(*fptr)(double, double))){
  MatrixXd c1 = MatrixXd::Zero(c0.rows(), c0.cols());
  MatrixXd I = MatrixXd::Identity(c0.rows(), c0.cols());
  VectorXd b0, b1, s, y, p, x1;
  double lambda = 0;
  double r = 0;
  int k;

  //Schritt 1
  b0 = gradient(x0, f);

  //Schritt 2: Liniensuchschritt
  lambda = newton(x0, b0, lambda, f);
  x1 = x0 + lambda * b0;
  b1 = gradient(x1, f);

  s = x1 - x0;
  y = b1 - b0;

  //Schritt 3: bereits gegeben durch c0

  k = 0;

  //Schritt 4 und 5
  do{
    //Schritt 4
    r = 1/(y.dot(s));
    c1 = (I - r * s * y.transpose()) * c0 * (I - r * y * s.transpose()) + r * s * s.transpose();
    k++;

    //Schritt 5
    p = -c1 * b1;
    x1 = x0 + p;
    b1 = gradient(x1, f);

    s = x1 - x0;
    y = b1 - b0;

    b0 = b1;
    c0 = c1;
    x0 = x1;

  } while(b1.norm() > epsilon);
  cout << "Iterationszahl: " << k << endl;
  cout << "f: " << f(x1(0), x1(1)) << endl;
  return x1;
}


int main(){
  VectorXd start(2);
    start << -1, 1;
  const double epsilon = 1e-5;
  cout << "Berechnung mit exakter Hessematrix" << endl;
  cout << BFGS(start, H_exakt(start, rosenbrock), epsilon, rosenbrock, gradient) << endl;
  cout << "Berechnung mit Diagonalmatrix" << endl;
  cout << BFGS(start, H_diagonal(start, rosenbrock), epsilon, rosenbrock, gradient) << endl;
  cout << "Berechnung mit Einheitsmatrix" << endl;
  cout << BFGS(start, H_einheitsmatrix(start, rosenbrock), epsilon, rosenbrock, gradient) << endl;
}
