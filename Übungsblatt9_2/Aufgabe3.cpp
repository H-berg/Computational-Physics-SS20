#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <random>
#include "Eigen/Dense"
#define _USE_MATH_DEFINES

using namespace Eigen;
using namespace std;

//Benennung von Dateien
string name(const string &basename, double d, double S, const string &ext)
{
	ostringstream name;
	name << basename << "_d_" << d << "_S_" << S << ext;
	return name.str();
}

//Berechne Strecke der Ortsvektoren für Permutation pi
double Strecke(MatrixXd &r, VectorXi &pi, const int N)
{
    double s = (r.col(pi(0)) - r.col(pi(N - 1))).norm();
    for (int i = 1; i < N; i++)
    {
        s += (r.col(pi(i)) - r.col(pi(i - 1))).norm();
    }
    return s;
}

//Lösung des Traveling Salesman Problems mit dem Metropolis Algorithmus
double TSP(MatrixXd &r, VectorXi &pi, double T_start, double T_end, double d, const int S, const int N)
{
    //Initialisierung
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<int> distribution_int(1, N - 2);//Der erste und letzte Punkt sollen festgehalten werden
    uniform_real_distribution<double> distribution_double(0., 1.);
    double s = Strecke(r, pi, N), s_vorschlag, s_min = Strecke(r, pi, N), diff, p, T = T_start; // k_B=1
    int temp, a, b;
    VectorXi pi_opt = pi, pi_vorschlag;

    while (T >= T_end)
    {
        for (int i = 0; i <= S; i++)
        {
            pi_vorschlag = pi;

            //Tausche zwei beliebige Orte a und b -> initialisiere Vorschlags Permutationsvektor
            a = distribution_int(generator);
            b = distribution_int(generator);
            //a und b sollen nicht gleich sein
            while (a == b)
            {
                b = distribution_int(generator);
            }
            pi_vorschlag.row(a).swap(pi_vorschlag.row(b));

            //Berechne erneut Strecke
            s_vorschlag = Strecke(r, pi_vorschlag, N);

            //Entscheide über Annahme oder Ablehung
            diff = s_vorschlag - s;
            if (diff <= 0)
            {
                pi = pi_vorschlag;
                s = s_vorschlag;

                //Ist das die bisher kürzeste Strecke?
                if (s < s_min)
                {
                    s_min = s;
                    pi_opt = pi;
                }
            }
            else
            {
                p = distribution_double(generator);
                if (p < exp(-diff/T))
                {
                    pi = pi_vorschlag;
                    s = s_vorschlag;
                }
            }
        }
        T *= d;
    }
    //Speicher die Konfiguration

    ofstream f;
    f.open(name("Strecke", d, S, ".txt"));
    f << "#x \t y \t Strecke: " << s_min << "\n";
    for (int i = 0; i < N; i++)
    {
        f << r(0, pi_opt(i)) << "\t" << r(1, pi_opt(i)) << "\n";
    }

    f.flush();
    f.close();
    
    //pi_opt soll zum schluss in pi stehen und die kürzeste Strecke zurückgegeben werden
    pi = pi_opt;
    return s_min;
}


//Initialisierung der Startvektoren gemäß Abbildung 1
void Initialisierung(MatrixXd &r, VectorXi &pi, const double Delta, const int N)
{
    VectorXd r_init = VectorXd::Zero(2);
    r.col(0) = r_init;
    double count = 1;

    //linke Kante
    while (r_init(1) < 3)
    {
        r_init(1) += Delta;
        r.col(count) = r_init;
        count++;
    }
    //Oben
    while (r_init(0) < 4)
    {
        r_init(0) += Delta;
        r.col(count) = r_init;
        count++;
    }
    //rechte Kante
    while (r_init(1) >= 0.2)
    {
        r_init(1) -= Delta;
        r.col(count) = r_init;
        count++;
    }
    //Korrektur*
    r_init(1) = 0;
    r(1, count-1)=0;
    //unten rechts
    while (r_init(0) > 3)
    {
        r_init(0) -= Delta;
        r.col(count) = r_init;
        count++;
    }
    //Innen rechts
    while (r_init(1) <= 1.8)
    {
        r_init(1) += Delta;
        r.col(count) = r_init;
        count++;
    }
    //Innen Oben
    while (r_init(0) > 1)
    {
        r_init(0) -= Delta;
        r.col(count) = r_init;
        count++;
    }
    //Innen links
    while (r_init(1) >= 0.2)
    {
        r_init(1) -= Delta;
        r.col(count) = r_init;
        count++;
    }
    //Korrektur*
    r_init(1)=0;
    r(1, count-1)=0;
    //Innen Oben
    while (r_init(0) >= 0.2)
    {
        r_init(0) -= Delta;
        r.col(count) = r_init;
        count++;
    }
    //*Aufgrund der vielen Subtraktionen von Delta kommt es zu kleinen numerischen Unsicherheiten, sodass r nicht ganz auf null zurück geht,
    //sondern nur auf 1E-17, das wurde nach den Schleifen die jeweils bis 0 laufen korrigiert

    
    //Ortsvektoren in der Optimalen Konfiguration abspeichern
    ofstream f;
    f.open("Initialisierung.txt");
    f << "x \t y \n";
    for (int i = 0; i < N; i++)
    {
        f << r(0, i) << "\t" << r(1, i) << "\n";
    }

    f.flush();
    f.close();

    //Vertausche Reihenfolge der Ortsvektoren
    pi = VectorXi::LinSpaced(N, 0, N-1);
    shuffle(pi.data()+1, pi.data()+pi.size()-1, default_random_engine(1908)); //Starts und Endvektor werden festgehalten

    //Speicher die
    f.open("Startkonfig.txt");
    f << "#x \t y \t s= " << Strecke(r, pi, N) <<  "\n";
    for (int i = 0; i < N; i++)
    {
        f << r(0, pi(i)) << "\t" << r(1, pi(i)) << "\n";
    }

    f.flush();
    f.close();
}

int main()
{
    const int N = 90;
    const double Delta = 0.2, T_start = 10, T_end = 1E-2;
    MatrixXd r = MatrixXd::Zero(2, N);
    VectorXi pi = VectorXi::Zero(N);
    double s;

    Initialisierung(r, pi, Delta, N);

    //Für verschiedene Werte für d und S wird die kürzeste Strecke berechnet
    //d = 0.9, S = 10
    s = TSP(r, pi, T_start, T_end, 0.9, 10, N);
     cout << "d = 0.9 \t S = 10 \t Kürzeste Strecke: " << s << endl;

    //d = 0.9, S = 100
    s = TSP(r, pi, T_start, T_end, 0.9, 100, N);
    cout << "d = 0.9 \t S = 100 \t Kürzeste Strecke: " << s << endl;

    //d = 0.9, S = 1000
    s = TSP(r, pi, T_start, T_end, 0.9, 1000, N);
    cout << "d = 0.9 \t S = 1000 \t Kürzeste Strecke: " << s << endl;

    //d = 0.9, S = 10000
    s = TSP(r, pi, T_start, T_end, 0.9, 10000, N);
    cout << "d = 0.9 \t S = 10000 \t Kürzeste Strecke: " << s << endl;

    //d = 0.99, S = 10
    s = TSP(r, pi, T_start, T_end, 0.99, 10, N);
    cout << "d = 0.99 \t S = 10 \t Kürzeste Strecke: " << s << endl;

    //d = 0.99, S = 100
    s = TSP(r, pi, T_start, T_end, 0.99, 100, N);
    cout << "d = 0.99 \t S = 100 \t Kürzeste Strecke: " << s << endl;

    //d = 0.99, S = 1000
    s = TSP(r, pi, T_start, T_end, 0.99, 1000, N);
    cout << "d = 0.99 \t S = 1000 \t Kürzeste Strecke: " << s << endl;
    
    /* läuft sehr lange wurde daher auskommentiert, lief aber bereits einmal durch
    //d = 0.99, S = 10000
    s = TSP(r, pi, T_start, T_end, 0.99, 10000, N);
    cout << "d = 0.99 \t S = 10000 \t Kürzeste Strecke: " << s << endl;
    */

    //d = 0.999, S = 10
    s = TSP(r, pi, T_start, T_end, 0.999, 10, N);
    cout << "d = 0.999 \t S = 10 \t Kürzeste Strecke: " << s << endl;

    //d = 0.999, S = 100
    s = TSP(r, pi, T_start, T_end, 0.999, 100, N);
    cout << "d = 0.999 \t S = 100 \t Kürzeste Strecke: " << s << endl;

    //d = 0.999, S = 1000
    s = TSP(r, pi, T_start, T_end, 0.999, 1000, N);
    cout << "d = 0.999 \t S = 1000 \t Kürzeste Strecke: " << s << endl;

    /* läuft sehr lange wurde daher auskommentiert, lief aber bereits einmal durch
    //d = 0.999, S = 10000
    s = TSP(r, pi, T_start, T_end, 0.999, 10000, N);
    cout << "d = 0.999 \t S = 10000 \t Kürzeste Strecke: " << s << endl;
    */
    return 0;
}