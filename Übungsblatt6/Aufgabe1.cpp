#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include "Eigen/Dense"
#define _USE_MATH_DEFINES

//using namespace std;
using namespace Eigen;
using namespace std;

//numerische Berechnung der Ableitung
//Die Ableitung wird durch den Differentialquotienten genähert
double diff(double (*fkt)(VectorXd), VectorXd x, VectorXd g, double lambda)
{
    double h = 1E-3;
    double r = (fkt(x + (lambda + h) * g) - fkt(x + (lambda - h) * g)) / (2 * h);
    return r;
}

double diffdiff(double (*fkt)(VectorXd), VectorXd x, VectorXd g, double lambda)
{
    double h = 1E-3;
    double r = (fkt(x + (lambda + h) * g) - 2 * fkt(x + lambda * g) + fkt(x + (lambda - h) * g)) / (h * h);

    return r;
}

//Newton verfahren lambda wird als pointer übergeben, damit dieser gleich dem Lambda entspricht welches die Fkt minimiert
void Newton(double (*fkt)(VectorXd), VectorXd x, double x_c, double &lambda, VectorXd g)
{
    double deltax;
    do
    {
        deltax = diff(fkt, x, g, lambda) / diffdiff(fkt, x, g, lambda);
        lambda -= deltax;
    } while (deltax > x_c);
}

VectorXd Gradientenverfahren(VectorXd x_0, double (*fkt)(VectorXd), VectorXd (*grad)(VectorXd), vector<double> &X) //ich verwende hier einen Std Vektor da es den Befehl append gibt
{
    VectorXd x = x_0;
    VectorXd g;
    double lambda;

    do
    {
        X.push_back(x(0));
        X.push_back(x(1));
        g = -grad(x);
        Newton(fkt, x, 1E-9, lambda, g);
        x = x + (lambda * g);
    } while (g.norm() > 1E-3);

    return x;
}

VectorXd konjugiert(VectorXd x_0, double (*fkt)(VectorXd), VectorXd (*grad)(VectorXd), vector<double> &X)
{
    VectorXd x = x_0;
    VectorXd g = -grad(x);
    VectorXd g_temp;
    VectorXd p = g;
    double mu;
    double lambda;

    do
    {
        lambda = 0; // lambda ist der Startwert des newtonverfahrens und muss daher in jeder Iteration erneut auf 0 gesetzt werden

        X.push_back(x(0));
        X.push_back(x(1));
        Newton(fkt, x, 1E-9, lambda, p);
        x = x + (lambda * p);
        g_temp = g;
        g = -grad(x);
        mu = g.squaredNorm() / (g_temp.squaredNorm());
        p = g + mu * p;
    } while (g.norm() > 1E-3);

    return x;
} 

//Funktion zum Aufgabenteil a)
double Fkt_1(VectorXd x)
{
    double r = (1 - x(0)) * (1 - x(0)) + 100 * (x(1) - x(0) * x(0)) * (x(1) - x(0) * x(0));
    return r;
}

//Gradient zur Funktion in Aufgabenteil a)
VectorXd DiffFkt_1(VectorXd x)
{
    VectorXd r(2);
    r(0) = -2 * (1 - x(0)) - 400 * x(0) * (x(1) - x(0) * x(0));
    r(1) = 200 * (x(1) - x(0) * x(0));
    return r;
}

double Fkt_2(VectorXd x)
{
    double r = 1 / (1 + ((exp(-10 * (x(0) * x(1) - 3) * (x(0) * x(1) - 3))) / (x(0) * x(0) + x(1) * x(1))));
    return r;
}

VectorXd DiffFkt_2(VectorXd x)
{
    VectorXd r(2);
    double e = exp(10 * (x(0) * x(1) - 3) * (x(0) * x(1) - 3))+1;
    double a = (e * (x(0) * x(0) + x(1) * x(1))) * (e * (x(0) * x(0) + x(1) * x(1)));
    r(0) = (2 * e * (10 * pow(x(0), 3) * pow(x(1), 2) - 30 * pow(x(0), 2) * x(1) + 10 * x(0) * pow(x(1), 4) + x(0) - 30 * pow(x(1), 3))) / a;
    r(1) = (2 * e * (10 * pow(x(0), 4) * x(1) - 30 * pow(x(0), 3) + 10 * x(0) * x(0) * pow(x(1), 3) - 30 * pow(x(1), 2) * x(0) + x(1))) / a;
    return r;
}

int main()
{
    VectorXd x_0(2);
    x_0 << -1, 1;
    vector<double> X, X2;

    //Die berechneten Minima werden zur Kontrolle ausgegeben, die einzelnen Schritte befinden sich danach im Vektor X
    cout << Gradientenverfahren(x_0, Fkt_1, DiffFkt_1, X) << endl;
    cout << konjugiert(x_0, Fkt_1, DiffFkt_1, X2) << endl;

    //Fehler
    double err;

    ofstream f;
    f.open("data_1.txt");
    f << "#Gradientenverfahren \n ";
    for (int i = 0; i < X.size(); i = i + 2)
    {
        err = sqrt((X[i] - 1) * (X[i] - 1) + (X[i + 1] - 1) * (X[i + 1] - 1)); // L_2 Norm zwischen dem Minimum (1, 1) und dem Vektor
        f << X[i] << "\t" << X[i + 1] << "\t" << err << "\n";
    }
    f.flush();
    f.close();

    f.open("data_2.txt");
    f << "#konjugierten Gradientenverfahren \n";
    for (int i = 0; i < X2.size(); i = i + 2)
    {
        err = sqrt((X2[i] - 1) * (X2[i] - 1) + (X2[i + 1] - 1) * (X2[i + 1] - 1));
        f << X2[i] << "\t" << X2[i + 1] << "\t" << err << "\n";
    }
    f.flush();
    f.close();

    //Aufgabenteil b)

    VectorXd x_1(2);
    VectorXd x_2(2);
    VectorXd x_3(2);
    x_1 << 1.5, 2.5;
    x_2 << -1.7, -1.9;
    x_3 << 0.5, 0.6;
    vector<double> X3;

    //Die berechneten Minima werden zur Kontrolle ausgegeben, die einzelnen Schritte befinden sich danach im Vektor X
    cout << konjugiert(x_1, Fkt_2, DiffFkt_2, X3) << endl;
    cout << konjugiert(x_2, Fkt_2, DiffFkt_2, X3) << endl;
    cout << konjugiert(x_3, Fkt_2, DiffFkt_2, X3) << endl;


    return 0;
}